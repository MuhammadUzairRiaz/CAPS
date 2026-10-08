using System.Globalization;
using System.Text.RegularExpressions;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>What the render overlays draw on an image of Width × Height: a text label (lines filled from live
/// attributes), a vertical colour legend, a true scale bar and an axis tripod.</summary>
public sealed record RenderSpec(string[] Lines, bool LabelBox, bool Legend, string LegendLo, string LegendHi, string LegendName,
                                double BarAngstrom, double BarPx, bool Tripod, double Yaw, double Pitch, uint Ink, bool Dark, double Width, double Height,
                                System.Text.Json.Nodes.JsonArray? Custom = null,   // a Python overlay's drawing (caps.overlay commands)
                                (string Label, uint Rgb)[]? LegendEntries = null)   // a categorical legend (the pipeline's colour coding)
{
    public double Unit => Width / 900.0;                         // the board draws a 1920-wide frame about 900 px across
    public string BarLabel => BarAngstrom.ToString("0.##", CultureInfo.InvariantCulture) + " Å";
}

/// <summary>Studio › Render (design/boards/RenderOverlays): the view with a render-frame guide, overlays that fill from
/// the frame's attributes, ambient occlusion, and an image or an image sequence written at the output size.</summary>
public sealed partial class MainViewModel
{
    public bool IsRender => _module == 19;
    public bool ShowAnalysisPanel => _module is not (19 or 20);
    public bool IsAnalyzeRail => _module is 1 or 20 or 22 or 23 or 24 or 28;

    public static readonly string[] RenderAntialias = ["1 sample", "4× (2 × 2)", "9× (3 × 3)", "16× (4 × 4)"];
    public static readonly string[] RenderFrameChoices = ["This frame", "All frames", "Every 10th frame"];
    public static readonly string RenderTokens = "[Title] [SourceFrame] [Frames] [Particles] [Molecules] [Density] [CellA] [CellB] [CellC] [Volume] [Mass]";

    private decimal _renderW = 1920, _renderH = 1080;
    private int _renderBg, _renderAa = 1, _renderFrames;
    private bool _renderAo = true, _renderDepth = true, _ovLabel = true, _ovLegend = true, _ovBar = true, _ovTripod = true;
    private string _ovText = "[Title] · frame [SourceFrame] · [Particles] atoms", _ovSub = "ρ = [Density] g/cm³ · L = [CellA] Å";
    private bool _rendering;
    private volatile bool _renderStop;
    private string _renderProgress = "";

    public void OpenRender()
    {
        if (_doc == null) { Status = "Open or build a structure first"; return; }
        SetModule(19);
        RenderRequested?.Invoke();
    }

    public decimal RenderW { get => _renderW; set { if (Set(ref _renderW, Math.Clamp(Math.Round(value), 64, 8192))) RenderChanged(); } }
    public decimal RenderH { get => _renderH; set { if (Set(ref _renderH, Math.Clamp(Math.Round(value), 64, 8192))) RenderChanged(); } }
    public int RenderBg { get => _renderBg; set { if (Set(ref _renderBg, value)) { Raise(nameof(RenderBgDark)); Raise(nameof(RenderBgWhite)); Raise(nameof(RenderBgTransparent)); RenderChanged(); } } }
    public bool RenderBgDark { get => _renderBg == 0; set { if (value) RenderBg = 0; } }
    public bool RenderBgWhite { get => _renderBg == 1; set { if (value) RenderBg = 1; } }
    public bool RenderBgTransparent { get => _renderBg == 2; set { if (value) RenderBg = 2; } }
    public bool RenderAo { get => _renderAo; set { if (Set(ref _renderAo, value)) RenderChanged(); } }
    public bool RenderDepth { get => _renderDepth; set { if (Set(ref _renderDepth, value)) RenderChanged(); } }
    public int RenderAa { get => _renderAa; set => Set(ref _renderAa, Math.Clamp(value, 0, 3)); }
    public int RenderFrames { get => _renderFrames; set => Set(ref _renderFrames, value); }
    public bool OvLabel { get => _ovLabel; set { if (Set(ref _ovLabel, value)) RenderChanged(); } }
    public string OvText { get => _ovText; set { if (Set(ref _ovText, value)) RenderChanged(); } }
    public string OvSub { get => _ovSub; set { if (Set(ref _ovSub, value)) RenderChanged(); } }
    public bool OvLegend { get => _ovLegend; set { if (Set(ref _ovLegend, value)) RenderChanged(); } }
    public bool OvBar { get => _ovBar; set { if (Set(ref _ovBar, value)) RenderChanged(); } }
    public bool OvTripod { get => _ovTripod; set { if (Set(ref _ovTripod, value)) RenderChanged(); } }
    public string OvLegendNote => _pipeHasLegend && PipelineRows.Count > 0 ? $"the Visualize colour coding · {_pipeLegendTitle}"
        : _colour == 3 ? "distance to molecule centre · viridis" : "shown when colouring by a property or by a Visualize colour coding";
    public string OvBarNote => _perspective ? "true only at the cell centre in perspective" : "true length · orthographic";
    public string OvLabelNote => ResolveTokens(_ovText, _frame);
    public bool Rendering { get => _rendering; private set { if (Set(ref _rendering, value)) Raise(nameof(RenderIdle)); } }
    public bool RenderIdle => !_rendering;
    public string RenderProgress { get => _renderProgress; private set => Set(ref _renderProgress, value); }
    public string RenderSizeText => $"render frame · {_renderW:0} × {_renderH:0}";

    public event Action? RenderOverlayChanged;
    private void RenderChanged()
    {
        Raise(nameof(RenderSizeText));
        Raise(nameof(OvLabelNote));
        RenderRequested?.Invoke();
        RenderOverlayChanged?.Invoke();
    }

    /// <summary>[Token]s filled from the frame: title, frame number, atom and molecule counts, density, cell edges.</summary>
    public string ResolveTokens(string template, int frame)
    {
        if (_doc == null || template.Length == 0) return template;
        var s = _doc.Summary();
        var inv = CultureInfo.InvariantCulture;
        return Regex.Replace(template, @"\[(\w+)\]", m => m.Groups[1].Value switch
        {
            "Title" => Title.Replace(" (unsaved)", ""),
            "SourceFrame" or "Frame" => frame.ToString(inv),
            "Frames" => s.Frames.ToString(inv),
            "Particles" or "Atoms" => s.Atoms.ToString("N0", inv),
            "Molecules" => s.Molecules.ToString("N0", inv),
            "Density" => s.CellValid != 0 ? s.Density.ToString("F3", inv) : "–",
            "CellA" => s.CellValid != 0 ? s.CellA.ToString("F1", inv) : "–",
            "CellB" => s.CellValid != 0 ? s.CellB.ToString("F1", inv) : "–",
            "CellC" => s.CellValid != 0 ? s.CellC.ToString("F1", inv) : "–",
            "Volume" => s.CellValid != 0 ? s.Volume.ToString("N0", inv) : "–",
            "Mass" => s.TotalMass.ToString("N1", inv),
            _ => m.Value,
        });
    }

    /// <summary>Render options of the output image (w × h) for background bg.</summary>
    public CapsRenderOpts RenderOptionsFor(int w, int h, int supersample)
    {
        var o = ViewOptions(w, h, supersample);
        o.LodNear = o.LodFar = 0;   // renders in full detail
        o.Background = _renderBg;
        o.Highlight0 = o.Highlight1 = o.Highlight2 = o.Highlight3 = -1;
        o.Focus = 0;
        o.AmbientOcclusion = _renderAo ? 1 : 0;
        o.DepthCue = _renderDepth ? 1 : 0;
        return o;
    }

    /// <summary>The overlay for an image of w × h pixels at pxPerAngstrom; dark marks whether the background is dark.</summary>
    public RenderSpec RenderSpecFor(int w, int h, double pxPerAngstrom, bool dark, int frame, System.Text.Json.Nodes.JsonArray? custom = null)
    {
        var lines = _ovLabel ? new[] { ResolveTokens(_ovText, frame), ResolveTokens(_ovSub, frame) }.Where(l => l.Length > 0).ToArray() : [];
        double barA = 0, barPx = 0;
        if (_ovBar && pxPerAngstrom > 0)
        {
            var target = 0.1 * w / pxPerAngstrom;
            barA = new[] { 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000.0 }.LastOrDefault(x => x <= target);
            if (barA <= 0) barA = 1;
            barPx = barA * pxPerAngstrom;
        }
        // the legend: the Visualize pipeline's colour coding when it has one (its property and range, or its categories),
        // else the view's distance-to-centre colouring
        if (_pipeHasLegend && PipelineRows.Count > 0)
        {
            (string, uint)[]? cats = null;
            if (!_pipeLegendContinuous)
                cats = PipeLegendEntries.Select(e => (e.Label, e.Brush is Avalonia.Media.ISolidColorBrush b ? (uint)(b.Color.ToUInt32() & 0xFFFFFF) : 0x888888u)).Take(12).ToArray();
            return new RenderSpec(lines, dark, _ovLegend, _pipeLegendLo, _pipeLegendHi, _pipeLegendTitle, barA, barPx, _ovTripod,
                                  Camera.Yaw, Camera.Pitch, dark ? 0xE9ECEFu : 0x141413u, dark, w, h, custom ?? OverlayForGuide(frame, w, h, dark), cats);
        }
        return new RenderSpec(lines, dark, _ovLegend && _colour == 3, LegendLo, LegendHi, "DistanceToCOM", barA, barPx, _ovTripod,
                              Camera.Yaw, Camera.Pitch, dark ? 0xE9ECEFu : 0x141413u, dark, w, h, custom ?? OverlayForGuide(frame, w, h, dark));
    }

    // ---------------------------------------------------------------- Python overlay (caps.overlay)
    // A script draws on the render with a small canvas (text, lines, shapes, inset plots) from the frame's attributes and
    // tables; python3 -m caps.overlay runs it and returns the drawing as commands, painted over the image and the guide.
    private bool _ovPython;
    private string _ovScript = "", _ovPythonNote = "", _ovConsole = "";
    private readonly Dictionary<(int Frame, int W, int H, bool Dark), System.Text.Json.Nodes.JsonArray> _ovCache = new();
    private readonly HashSet<(int, int, int, bool)> _ovPending = new();
    private string _ovCacheKey = "";
    public bool OvPython { get => _ovPython; set { if (Set(ref _ovPython, value)) { _ovCache.Clear(); RenderChanged(); } } }
    public string OvScript { get => _ovScript; set { if (Set(ref _ovScript, value)) { _ovCache.Clear(); Raise(nameof(OvScriptName)); WatchOverlay(); RenderChanged(); } } }
    private FileSystemWatcher? _ovWatch;
    /// <summary>Saving the script redraws the guide (the cache is keyed by its modification time).</summary>
    private void WatchOverlay()
    {
        _ovWatch?.Dispose();
        _ovWatch = null;
        if (_ovScript.Length == 0 || System.IO.Path.GetDirectoryName(_ovScript) is not { Length: > 0 } dir || !Directory.Exists(dir)) return;
        try
        {
            _ovWatch = new FileSystemWatcher(dir, System.IO.Path.GetFileName(_ovScript)) { NotifyFilter = NotifyFilters.LastWrite | NotifyFilters.Size, EnableRaisingEvents = true };
            _ovWatch.Changed += (_, _) => Avalonia.Threading.Dispatcher.UIThread.Post(() => RenderOverlayChanged?.Invoke());
        }
        catch { /* no watching on this file system: the next change of the view redraws */ }
    }
    public string OvScriptName => _ovScript.Length == 0 ? "no script yet" : System.IO.Path.GetFileName(_ovScript);
    public string OvPythonNote { get => _ovPythonNote; private set => Set(ref _ovPythonNote, value); }
    public string OvConsole { get => _ovConsole; private set => Set(ref _ovConsole, value); }
    public static string OverlayFolder => AppSettings.Override != null ? System.IO.Path.Combine(System.IO.Path.GetDirectoryName(AppSettings.Override)!, "caps-overlays")
                                                                      : System.IO.Path.Combine(AppSettings.Folder, "overlays");

    private const string OverlayTemplate = """
        # A CAPS render overlay: draw on the image with the canvas (pixels, y down from the top left).
        # data.attributes: the Visualize pipeline's global attributes on this frame; data.tables: its data tables.
        from caps.overlay import overlay


        @overlay
        def draw(canvas, data):
            rho = data.attributes.get("Density", 0)
            canvas.text(24, canvas.height - 64, f"frame {data.frame} · ρ = {rho:.3f} g/cm³", size=26)
            # an inset plot of a pipeline table, e.g. the g(r) of Coordination & RDF
            rdf = data.tables.get("rdf")
            if rdf:
                canvas.plot(rdf.column(0), rdf.column(1), box=(canvas.width - 520, 40, 480, 300), title="g(r)", xlabel="r (Å)")

        """;

    /// <summary>A new overlay script from the template in ~/.caps/overlays (kept when it exists), chosen.</summary>
    public string NewOverlayScript()
    {
        Directory.CreateDirectory(OverlayFolder);
        var path = System.IO.Path.Combine(OverlayFolder, "overlay.py");
        if (!File.Exists(path)) File.WriteAllText(path, OverlayTemplate);
        OvScript = path;
        OvPython = true;
        return path;
    }

    /// <summary>What the overlay script is given for a frame: size, frame, title, the pipeline's attributes and tables.</summary>
    private static string OverlayInput(CapsDocument doc, int frame, int frames, string title, int w, int h, bool dark)
    {
        var o = new System.Text.Json.Nodes.JsonObject { ["width"] = w, ["height"] = h, ["dark"] = dark, ["frame"] = frame, ["frames"] = frames, ["title"] = title };
        var attrs = new System.Text.Json.Nodes.JsonObject();
        var tables = new System.Text.Json.Nodes.JsonObject();
        string result;
        try { result = doc.PipelineResult(); } catch { result = ""; }
        if (result.Length > 0 && System.Text.Json.Nodes.JsonNode.Parse(result) is System.Text.Json.Nodes.JsonObject r)
        {
            if (r["attributes"] is System.Text.Json.Nodes.JsonArray a)
                foreach (var x in a) if ((string?)x?["name"] is { } n) attrs[n] = (double?)x!["value"] ?? 0;
            if (r["tables"] is System.Text.Json.Nodes.JsonArray ts)
                foreach (var t in ts.OfType<System.Text.Json.Nodes.JsonObject>())
                    if ((string?)t["name"] is { } n) tables[n] = new System.Text.Json.Nodes.JsonObject { ["columns"] = t["columns"]?.DeepClone(), ["rows"] = t["rows"]?.DeepClone() };
        }
        else
        {
            var s = doc.Summary();
            attrs["Particles"] = s.Atoms; attrs["Molecules"] = s.Molecules; attrs["Mass"] = s.TotalMass;
            if (s.CellValid != 0) { attrs["Density"] = s.Density; attrs["CellVolume"] = s.Volume; }
        }
        o["attributes"] = attrs;
        o["tables"] = tables;
        return o.ToJsonString();
    }

    /// <summary>Runs the overlay script on one frame's data: its commands, what it printed, and an error when it failed.</summary>
    private static (System.Text.Json.Nodes.JsonArray? Commands, string Console, string Error) RunOverlay(string script, string input)
    {
        try
        {
            var python = PythonExe;   // Settings › Python & scripting, else CAPS_PYTHON, else python3
            var psi = new System.Diagnostics.ProcessStartInfo(python)
            {
                RedirectStandardInput = true, RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false,
                WorkingDirectory = System.IO.Path.GetDirectoryName(script) ?? ".",
            };
            PythonProcess.Utf8Io(psi);
            psi.ArgumentList.Add("-m");
            psi.ArgumentList.Add("caps.overlay");
            psi.ArgumentList.Add(script);
            if (Paths.Python is { } pkg) psi.Environment["PYTHONPATH"] = pkg + (Environment.GetEnvironmentVariable("PYTHONPATH") is { Length: > 0 } pp ? System.IO.Path.PathSeparator + pp : "");
            psi.Environment["CAPS_LIB"] = NativeLibraryPath;
            psi.StandardInputEncoding = new System.Text.UTF8Encoding(false);
            psi.StandardOutputEncoding = System.Text.Encoding.UTF8;
            using var proc = System.Diagnostics.Process.Start(psi) ?? throw new InvalidOperationException("cannot start " + python);
            proc.StandardInput.Write(input);
            proc.StandardInput.Close();
            var err = proc.StandardError.ReadToEndAsync();
            var outText = proc.StandardOutput.ReadToEnd();
            if (!proc.WaitForExit(30000)) { try { proc.Kill(true); } catch { } return (null, "", "the overlay script took longer than 30 s"); }
            var console = err.Result.Trim();
            if (proc.ExitCode != 0) return (null, console, console.Split('\n').LastOrDefault(l => l.Trim().Length > 0) ?? $"exit code {proc.ExitCode}");
            var j = System.Text.Json.Nodes.JsonNode.Parse(outText);
            return (j?["commands"] as System.Text.Json.Nodes.JsonArray ?? [], console, "");
        }
        catch (Exception e) { return (null, "", e.Message + " (choose a Python 3 interpreter in Settings › Python & scripting)"); }
    }

    /// <summary>The overlay's drawing for the guide in the view: from the cache, else worked out in the background.</summary>
    private System.Text.Json.Nodes.JsonArray? OverlayForGuide(int frame, int w, int h, bool dark)
    {
        if (!_ovPython || _ovScript.Length == 0 || _doc == null) return null;
        var stamp = File.Exists(_ovScript) ? File.GetLastWriteTimeUtc(_ovScript).Ticks.ToString(CultureInfo.InvariantCulture) : "";
        var key = $"{_ovScript}|{stamp}|{_pipeGen}";
        if (key != _ovCacheKey) { _ovCache.Clear(); _ovPending.Clear(); _ovCacheKey = key; }
        var k = (frame, w, h, dark);
        if (_ovCache.TryGetValue(k, out var hit)) return hit;
        if (!_ovPending.Add(k)) return null;
        var doc = _doc;
        var (script, frames, title) = (_ovScript, _frames, Title);
        var input = OverlayInput(doc, frame, frames, title, w, h, dark);
        Task.Run(() => RunOverlay(script, input)).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            _ovPending.Remove(k);
            var (cmds, console, error) = t.Result;
            OvConsole = console;
            OvPythonNote = error.Length > 0 ? "Overlay: " + error : $"{cmds?.Count ?? 0} drawing commands · frame {frame}";
            if (cmds != null && key == _ovCacheKey) _ovCache[k] = cmds;
            RenderOverlayChanged?.Invoke();
        }));
        return null;
    }

    /// <summary>The camera of the 3D view: in Render the view zooms so the output frame fills 86 % of it.</summary>
    public CapsCamera ViewCamera(double vw, double vh)
    {
        if (!IsRender || _doc == null || vw < 20 || vh < 20) return Camera;
        var (w, h) = ((int)_renderW, (int)_renderH);
        try
        {
            var sView = _doc.ViewScale(Camera, ViewOptions((int)vw, (int)vh, 1));
            var sOut = _doc.ViewScale(Camera, RenderOptionsFor(w, h, 1));
            if (sView <= 0 || sOut <= 0) return Camera;
            var k = sView / sOut;
            var f = Math.Min(0.86 * vw / (w * k), 0.86 * vh / (h * k));
            var c = Camera;
            c.Zoom *= f;
            // pan is in Å, unchanged by zoom
            return c;
        }
        catch { return Camera; }
    }

    /// <summary>Where the output falls in a vw × vh view: the frame guide (px) and the overlay drawn inside it.</summary>
    public (double X, double Y, double W, double H, RenderSpec Spec)? RenderGuide(double vw, double vh)
    {
        if (_doc == null || vw < 20 || vh < 20) return null;
        var (w, h) = ((int)_renderW, (int)_renderH);
        double sView, sOut;
        try
        {
            sView = _doc.ViewScale(ViewCamera(vw, vh), ViewOptions((int)vw, (int)vh, 1));
            sOut = _doc.ViewScale(Camera, RenderOptionsFor(w, h, 1));
        }
        catch { return null; }
        if (sView <= 0 || sOut <= 0) return null;
        var k = sView / sOut;   // view px per output px
        var fw = w * k;
        var fh = h * k;
        return ((vw - fw) / 2, (vh - fh) / 2, fw, fh, RenderSpecFor(w, h, sOut, _renderBg == 0 || _viewBackground == 0 && _renderBg == 2, _frame));
    }

    /// <summary>Renders frames (the current one, all, or every 10th) at the output size and hands each image to write
    /// (on the UI thread: it draws the overlay and saves). Returns the number written.</summary>
    public async Task<int> RenderOut(Func<int, string> pathOf, bool movie, Action<byte[], int, int, RenderSpec, string> write)
    {
        if (_doc == null || Rendering) return 0;
        var doc = _doc;
        var (w, h) = ((int)_renderW, (int)_renderH);
        var ss = _renderAa + 1;
        if ((long)w * h * ss * ss > 40_000_000) ss = Math.Max(1, (int)Math.Floor(Math.Sqrt(40_000_000.0 / w / h)));
        var frames = !movie || _renderFrames == 0 ? new List<int> { _frame }
                   : Enumerable.Range(0, _frames).Where(f => _renderFrames == 1 || f % 10 == 0).ToList();
        var cam = Camera;
        var opt = RenderOptionsFor(w, h, ss);
        var scale = doc.ViewScale(cam, opt);
        var dark = _renderBg == 0;
        Rendering = true;
        _renderStop = false;
        var keep = _frame;
        int done = 0;
        try
        {
            foreach (var f in frames)
            {
                if (_renderStop) break;
                RenderProgress = frames.Count > 1 ? $"frame {f} · {done + 1} of {frames.Count}" : $"rendering {w} × {h}";
                var rgba = new byte[w * h * 4];
                var (py, script, title, nf) = (_ovPython && _ovScript.Length > 0, _ovScript, Title, _frames);
                var drawn = await Task.Run(() =>
                {
                    if (frames.Count > 1) doc.SetFrame(f);
                    doc.Render(cam, opt, rgba);
                    return py ? RunOverlay(script, OverlayInput(doc, f, nf, title, w, h, dark)) : (null, "", "");
                });
                if (py)
                {
                    OvConsole = drawn.Console;
                    if (drawn.Error.Length > 0) OvPythonNote = "Overlay: " + drawn.Error;
                }
                write(rgba, w, h, RenderSpecFor(w, h, scale, dark, f, drawn.Commands), pathOf(f));
                done++;
            }
        }
        finally
        {
            if (frames.Count > 1 && _doc == doc) doc.SetFrame(keep);
            Rendering = false;
            RenderProgress = "";
            RenderRequested?.Invoke();
        }
        Status = done == 1 ? $"Rendered {w} × {h}{(_renderAo ? " · ambient occlusion" : "")} · {pathOf(frames[0])}"
                           : $"Rendered {done} frames at {w} × {h}{(_renderStop ? " (stopped)" : "")} · {System.IO.Path.GetDirectoryName(pathOf(frames[0]))}";
        return done;
    }

    public void StopRender() => _renderStop = true;
}
