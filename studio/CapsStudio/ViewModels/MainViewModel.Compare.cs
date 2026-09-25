using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using Avalonia.Media;
using Avalonia.Media.Imaging;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed record CompareRow(string Metric, string A, string B, string Diff, IBrush DiffBrush);

/// <summary>Analyze › Compare (design/boards/CompareCells): two inputs side by side with linked cameras, the Visualize
/// pipeline on both, every attribute as A, B and B − A, and a data table of both overlaid.</summary>
public sealed partial class MainViewModel
{
    public bool IsCompare => _module == 23;
    public ObservableCollection<CompareRow> CompareRows { get; } = new();
    public ObservableCollection<string> CompareTables { get; } = new();

    private CapsDocument? _cmpA, _cmpB;
    private string _cmpPathA = "", _cmpPathB = "", _cmpNote = "";
    private Bitmap? _cmpImgA, _cmpImgB;
    private bool _cmpLink = true;
    private int _cmpTable;
    private CapsCamera _camA = new() { Yaw = 0.55, Pitch = 0.40, Zoom = 1.0 }, _camB = new() { Yaw = 0.55, Pitch = 0.40, Zoom = 1.0 };
    private JsonNode? _cmpResA, _cmpResB;
    private int _cmpGen, _cmpReqA, _cmpReqB;
    public const int CompareTileW = 620, CompareTileH = 380;
    public event Action? ComparePlotChanged;

    public string ComparePathA => _cmpPathA;
    public string ComparePathB => _cmpPathB;
    public string CompareLabelA => "A · " + CompareName(_cmpPathA);
    public string CompareLabelB => "B · " + CompareName(_cmpPathB);
    public Bitmap? CompareImageA { get => _cmpImgA; private set => Set(ref _cmpImgA, value); }
    public Bitmap? CompareImageB { get => _cmpImgB; private set => Set(ref _cmpImgB, value); }
    public bool CompareHasB => _cmpB != null;
    public bool CompareLink { get => _cmpLink; set { if (Set(ref _cmpLink, value)) { if (value) _camB = _camA; RenderCompare(); } } }
    public string CompareNote { get => _cmpNote; private set => Set(ref _cmpNote, value); }
    public int CompareTable { get => _cmpTable; set { if (value >= 0 && Set(ref _cmpTable, value)) { Raise(nameof(CompareTableName)); ComparePlotChanged?.Invoke(); } } }
    /// <summary>The overlaid table by title (the picker binds to it, so it survives the list being rebuilt).</summary>
    public string? CompareTableName
    {
        get => _cmpTable < CompareTables.Count ? CompareTables[_cmpTable] : null;
        set { var i = value == null ? -1 : CompareTables.IndexOf(value); if (i >= 0) CompareTable = i; }
    }
    public string ComparePipelineText => BatchPipelineText;

    private static string CompareName(string p) => p.Length == 0 ? "choose a file" :
        $"{System.IO.Path.GetFileName(System.IO.Path.GetDirectoryName(p))}/{System.IO.Path.GetFileName(p)}";

    public void OpenCompare()
    {
        SetModule(23);
        Raise(nameof(ComparePipelineText));
        if (_cmpA == null && _doc?.Path is { Length: > 0 } p && File.Exists(p)) _ = SetCompareInput(true, p);
        else RefreshCompare();
    }

    /// <summary>Opens one side (off the UI thread; a dump picks up the data file beside it) and compares again.</summary>
    public async Task SetCompareInput(bool a, string path)
    {
        var topo = TopologyFor(path);
        var req = a ? ++_cmpReqA : ++_cmpReqB;
        try
        {
            var doc = await Task.Run(() => CapsDocument.Open(path, topo.Length > 0 ? topo : null));
            if (req != (a ? _cmpReqA : _cmpReqB)) { doc.Dispose(); return; }   // a later choice for this side won
            if (a) { _cmpA?.Dispose(); _cmpA = doc; _cmpPathA = path; }
            else { _cmpB?.Dispose(); _cmpB = doc; _cmpPathB = path; }
            foreach (var n in new[] { nameof(ComparePathA), nameof(ComparePathB), nameof(CompareLabelA), nameof(CompareLabelB), nameof(CompareHasB) }) Raise(n);
            RefreshCompare();
        }
        catch (Exception e) { CompareNote = $"Could not open {System.IO.Path.GetFileName(path)}: {e.Message}"; }
    }

    /// <summary>Runs the pipeline on both sides, fills the metric table and the common data tables, and renders both.</summary>
    public void RefreshCompare()
    {
        if (!IsCompare) return;
        var pipeline = PipelineRows.Count == 0 ? "" : PipelineJson();
        JsonNode? Run(CapsDocument? d)
        {
            if (d == null) return null;
            try { d.SetPipeline(pipeline); var r = d.PipelineResult(); return r.Length > 0 ? JsonNode.Parse(r) : null; }
            catch (Exception e) { CompareNote = e.Message; return null; }
        }
        _cmpResA = Run(_cmpA);
        _cmpResB = Run(_cmpB);
        var inv = CultureInfo.InvariantCulture;
        Dictionary<string, double> Attrs(JsonNode? r, CapsDocument? d)
        {
            var map = new Dictionary<string, double>();
            if (r?["attributes"] is JsonArray a) foreach (var x in a) map[(string?)x?["name"] ?? ""] = (double?)x?["value"] ?? double.NaN;
            else if (d != null) { var s = d.Summary(); map["Particles"] = s.Atoms; map["Bonds"] = s.Bonds; map["Molecules"] = s.Molecules; map["Density"] = s.Density; map["Mass"] = s.TotalMass; }
            return map;
        }
        var A = Attrs(_cmpResA, _cmpA);
        var B = Attrs(_cmpResB, _cmpB);
        CompareRows.Clear();
        var skip = new[] { "SourceFrame", "Timestep", "Selected" };
        foreach (var k in A.Keys.Concat(B.Keys).Distinct().Where(k => !skip.Contains(k)))
        {
            var hasA = A.TryGetValue(k, out var a) && double.IsFinite(a);
            var hasB = B.TryGetValue(k, out var b) && double.IsFinite(b);
            string F(double v) => Math.Abs(v) >= 1000 || Math.Abs(v - Math.Round(v)) < 1e-9 ? v.ToString("0.###", inv) : v.ToString("G4", inv);
            var diff = hasA && hasB ? b - a : double.NaN;
            var rel = hasA && hasB && Math.Abs(a) > 1e-12 ? Math.Abs(diff / a) : 0;
            CompareRows.Add(new CompareRow(k, hasA ? F(a) : "—", hasB ? F(b) : "—", double.IsFinite(diff) ? (diff >= 0 ? "+" : "") + F(diff) : "—",
                                           Tokens.Brush(rel > 0.1 ? "WarnB" : "MutedB")));
        }
        // data tables both sides have
        var titles = new List<string>();
        if (_cmpResA?["tables"] is JsonArray ta)
            foreach (var t in ta)
            {
                var title = (string?)t?["title"] ?? "";
                if (_cmpB == null || (_cmpResB?["tables"] is JsonArray tb && tb.Any(x => (string?)x?["title"] == title))) titles.Add(title);
            }
        if (!titles.SequenceEqual(CompareTables))
        {
            CompareTables.Clear();
            foreach (var t in titles) CompareTables.Add(t);
            _cmpTable = 0;
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { Raise(nameof(CompareTable)); Raise(nameof(CompareTableName)); });
        }
        CompareNote = _cmpB == null ? "Choose B to compare." :
            "Same pipeline on both. For two builds of one recipe with different seeds, the differences are seed-to-seed scatter: the error to expect from one build.";
        ComparePlotChanged?.Invoke();
        RenderCompare();
    }

    /// <summary>The chosen table's first two columns, for A and for B.</summary>
    public ((double, double)[] A, (double, double)[] B, string X, string Y) ComparePlot()
    {
        (double, double)[] Pick(JsonNode? r, out string x, out string y)
        {
            x = y = "";
            if (_cmpTable >= CompareTables.Count || r?["tables"] is not JsonArray ts) return [];
            var t = ts.FirstOrDefault(q => (string?)q?["title"] == CompareTables[_cmpTable]);
            if (t?["rows"] is not JsonArray rows || t["columns"] is not JsonArray cols || cols.Count < 2) return [];
            x = (string?)cols[0] ?? ""; y = (string?)cols[1] ?? "";
            return rows.Select(q => ((double?)q?[0] ?? 0, (double?)q?[1] ?? 0)).ToArray();
        }
        var a = Pick(_cmpResA, out var xl, out var yl);
        var b = Pick(_cmpResB, out _, out _);
        return (a, b, xl, yl);
    }

    public void RotateCompare(bool a, double dx, double dy)
    {
        ref var cam = ref a ? ref _camA : ref _camB;
        cam.Yaw += dx * 0.008;
        cam.Pitch = Math.Clamp(cam.Pitch + dy * 0.008, -Math.PI / 2, Math.PI / 2);
        if (_cmpLink) { if (a) _camB = _camA; else _camA = _camB; }
        RenderCompare();
    }

    private void RenderCompare()
    {
        var gen = ++_cmpGen;
        var opt = ViewOptions(CompareTileW, CompareTileH, 2);
        opt.Highlight0 = opt.Highlight1 = opt.Highlight2 = opt.Highlight3 = -1;
        opt.Focus = 0;
        foreach (var (doc, cam, isA) in new[] { (_cmpA, _camA, true), (_cmpB, _camB, false) })
        {
            if (doc == null) { if (isA) CompareImageA = null; else CompareImageB = null; continue; }
            var d = doc;
            var c = cam;
            Task.Run(() =>
            {
                var rgba = new byte[CompareTileW * CompareTileH * 4];
                try { d.Render(c, opt, rgba); } catch { return; }
                Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                {
                    if (gen != _cmpGen || !IsCompare) return;
                    var bmp = ToBitmap(rgba, CompareTileW, CompareTileH);
                    if (isA) CompareImageA = bmp; else CompareImageB = bmp;
                });
            });
        }
    }
}
