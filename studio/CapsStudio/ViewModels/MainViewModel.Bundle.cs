using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

public sealed record BundleRow(string File, string Content, string Size, string Hash);

/// <summary>Export › Figure bundle (design/boards/FigureBundle): the figure with the data behind it, the pipeline as run,
/// the input and a provenance record with every file's sha256; caps reproduce rebuilds it and compares.</summary>
public sealed partial class MainViewModel
{
    public bool IsBundle => _module == 26;
    public ObservableCollection<BundleRow> BundleRows { get; } = new();
    private bool _bInput, _bPipeline = true, _bData = true, _bReadme = true, _bBusy;
    private string _bName = "figure", _bProv = "", _bCsv = "", _bCsvName = "", _bNote = "";
    private int _bGen;

    public void OpenBundle()
    {
        if (_doc == null) { Status = "Open a structure first"; return; }
        if (_bName == "figure") _bName = System.IO.Path.GetFileNameWithoutExtension(_doc.Path is { Length: > 0 } p ? p : "figure") + "_figure";
        SetModule(26);
        Raise(nameof(BundleName)); Raise(nameof(BundleFileName)); Raise(nameof(BundleReproduce)); Raise(nameof(BundleInputAvailable));
        RefreshBundle();
    }

    public string BundleName { get => _bName; set { if (Set(ref _bName, value.Trim().Length == 0 ? "figure" : value.Trim())) { Raise(nameof(BundleFileName)); Raise(nameof(BundleReproduce)); RefreshBundle(); } } }
    public string BundleFileName => _bName + ".caps-bundle.zip";
    public bool BundleInput { get => _bInput; set { if (Set(ref _bInput, value)) RefreshBundle(); } }
    public bool BundlePipeline { get => _bPipeline; set { if (Set(ref _bPipeline, value)) RefreshBundle(); } }
    public bool BundleData { get => _bData; set { if (Set(ref _bData, value)) RefreshBundle(); } }
    public bool BundleReadme { get => _bReadme; set { if (Set(ref _bReadme, value)) RefreshBundle(); } }
    public bool BundleInputAvailable => _doc?.Path is { Length: > 0 } p && File.Exists(p);
    public string BundleProvenance { get => _bProv; private set => Set(ref _bProv, value); }
    public string BundleCsv { get => _bCsv; private set => Set(ref _bCsv, value); }
    public string BundleCsvName { get => _bCsvName; private set => Set(ref _bCsvName, value); }
    public string BundleNote { get => _bNote; private set => Set(ref _bNote, value); }
    public bool BundleBusy { get => _bBusy; private set { if (Set(ref _bBusy, value)) Raise(nameof(BundleIdle)); } }
    public bool BundleIdle => !_bBusy;
    public string BundleCount => $"{BundleRows.Count} files";
    public string BundleReproduce => $"caps reproduce {BundleFileName}";

    private string BundleOptions()
    {
        var input = _doc?.Path ?? "";
        var topo = input.Length > 0 ? TopologyFor(input) : "";
        var (w, h) = FigPixels;
        return new JsonObject
        {
            ["name"] = _bName, ["input"] = File.Exists(input) ? input : "", ["topology"] = topo, ["frame"] = _frame,
            ["width"] = w, ["height"] = h,
            ["include_input"] = _bInput && BundleInputAvailable, ["include_pipeline"] = _bPipeline, ["include_data"] = _bData, ["include_readme"] = _bReadme,
            ["pipeline"] = JsonNode.Parse(PipelineRows.Count == 0 ? "{\"steps\":[]}" : PipelineJson()),
        }.ToJsonString();
    }

    /// <summary>Lists the files the bundle will hold with the sizes and hashes known now (off the UI thread).</summary>
    public void RefreshBundle()
    {
        if (_doc == null || !IsBundle) return;
        var doc = _doc;
        var gen = ++_bGen;
        var opts = BundleOptions();
        BundleBusy = true;
        Task.Run(() =>
        {
            string? json = null, error = null;
            try { json = doc.BundlePreview(opts); } catch (Exception e) { error = e.Message; }
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                if (gen != _bGen) return;
                BundleBusy = false;
                BundleRows.Clear();
                if (json == null) { BundleNote = error ?? ""; return; }
                var j = JsonNode.Parse(json)!;
                var inv = CultureInfo.InvariantCulture;
                foreach (var f in (JsonArray)j["files"]!)
                {
                    var bytes = (double?)f!["bytes"] ?? -1;
                    var size = bytes < 0 ? "on export" : bytes < 1024 ? $"{bytes:0} B" : bytes < 1 << 20 ? (bytes / 1024).ToString("F1", inv) + " KB" : (bytes / 1048576).ToString("F1", inv) + " MB";
                    var hash = (string?)f["sha256"] ?? "";
                    BundleRows.Add(new BundleRow((string?)f["name"] ?? "", (string?)f["note"] ?? "", size, hash.Length >= 12 ? hash[..12] : "—"));
                }
                BundleProvenance = (string?)j["provenance"] ?? "";
                BundleCsv = (string?)j["csv"] ?? "";
                BundleCsvName = (string?)j["csv_name"] ?? "";
                BundleNote = PipelineRows.Count == 0 ? "The pipeline is empty: add steps in Visualize to put their data beside the figure." :
                    "Hashes are shown for files whose content is fixed now; the figures are made on export, with the camera of the view.";
                Raise(nameof(BundleCount));
            });
        });
    }

    public async Task<int> ExportBundle(string path)
    {
        if (_doc == null) return 0;
        var doc = _doc;
        var opts = BundleOptions();
        var cam = Camera;
        var (w, h) = FigPixels;
        var ro = ExportOptions(w, h);
        BundleBusy = true;
        try
        {
            var n = await Task.Run(() => doc.BundleWrite(path, opts, cam, ro));
            Status = $"Wrote {System.IO.Path.GetFileName(path)} · {n} files · reproduce with caps reproduce";
            return n;
        }
        catch (Exception e) { Status = "Bundle failed: " + e.Message; return 0; }
        finally { BundleBusy = false; }
    }
}
