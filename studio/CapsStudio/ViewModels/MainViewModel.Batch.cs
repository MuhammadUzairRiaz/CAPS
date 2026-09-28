using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json.Nodes;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One input of a batch: its state and the pipeline's attributes on it.</summary>
public sealed class BatchInput : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    public string Path { get; init; } = "";
    public string Topology { get; init; } = "";
    private string _shown = "";
    public string Shown { get => _shown.Length > 0 ? _shown : RecentFiles.Tilde(Path); set { _shown = value; Raise(nameof(Shown)); } }
    private string _state = "queued", _time = "—", _error = "";
    private string[] _values = [];
    public string State { get => _state; set { _state = value; Raise(nameof(State)); Raise(nameof(StateBrush)); } }
    public string Time { get => _time; set { _time = value; Raise(nameof(Time)); } }
    public string Error { get => _error; set { _error = value; Raise(nameof(Error)); } }
    public string[] Values { get => _values; set { _values = value; Raise(nameof(Values)); } }
    public Dictionary<string, double> Attributes { get; } = new();
    public IBrush StateBrush => Tokens.Brush(_state switch { "done" => "OkB", "failed" => "ErrB", "running" => "AccB", _ => "DimB" });
}

/// <summary>Analyze › Batch (design/boards/BatchRun): the Visualize pipeline on many inputs, one table out; failures are
/// skipped and kept as failed rows; results.csv records the pipeline's sha256 with every row.</summary>
public sealed partial class MainViewModel
{
    public bool IsBatch => _module == 22;
    public ObservableCollection<BatchInput> BatchInputs { get; } = new();
    public ObservableCollection<string> BatchColumns { get; } = new();
    public ObservableCollection<Row> BatchSummary { get; } = new();
    public static readonly string[] BatchFrames = ["Last frame", "First frame"];
    public static readonly string[] BatchWorkerChoices = ["1 local", "2 local", "4 local", "8 local"];

    private string _batchPattern = "", _batchState = "", _batchOut = "", _batchErrors = "";
    private int _batchFrame, _batchWorkers = 2, _batchMetric;
    private bool _batchRunning;
    private volatile bool _batchStop;
    private List<string> _batchKeys = new();

    public void OpenBatch()
    {
        SetModule(22);
        if (_batchPattern.Length == 0 && _doc?.Path is { Length: > 0 } p && File.Exists(p))
            BatchPattern = Path.Combine(Path.GetDirectoryName(Path.GetDirectoryName(p)) ?? "", "*", Path.GetFileName(p));
        Raise(nameof(BatchPipelineText));
    }

    public string BatchPattern { get => _batchPattern; set { if (Set(ref _batchPattern, value)) ExpandBatch(); } }
    public int BatchFrame { get => _batchFrame; set => Set(ref _batchFrame, value); }
    public int BatchWorkers { get => _batchWorkers; set => Set(ref _batchWorkers, Math.Clamp(value, 0, 3)); }
    public string BatchState { get => _batchState; private set => Set(ref _batchState, value); }
    public string BatchOut { get => _batchOut; private set => Set(ref _batchOut, value); }
    public string BatchErrors { get => _batchErrors; private set { if (Set(ref _batchErrors, value)) Raise(nameof(BatchHasErrors)); } }
    public bool BatchHasErrors => _batchErrors.Length > 0;
    public bool BatchRunning { get => _batchRunning; private set { if (Set(ref _batchRunning, value)) Raise(nameof(BatchIdle)); } }
    public bool BatchIdle => !_batchRunning;
    public int BatchMetric { get => _batchMetric; set { if (value >= 0 && Set(ref _batchMetric, value)) BatchPlotChanged?.Invoke(); } }
    public event Action? BatchPlotChanged;
    public string BatchPipelineText => PipelineRows.Count == 0 ? "No steps: the base attributes (density, particles, charge …)"
        : string.Join(" → ", PipelineRows.Reverse().Where(r => r.Enabled).Select(r => r.Title));

    /// <summary>The inputs matching the pattern: * and ? in any path segment (cells/seed*/PS_melt.data).</summary>
    private void ExpandBatch()
    {
        if (BatchRunning) return;
        BatchInputs.Clear();
        foreach (var f in Glob(_batchPattern).Take(500)) BatchInputs.Add(new BatchInput { Path = f, Topology = TopologyFor(f) });
        ShortenBatchPaths();
        BatchState = BatchInputs.Count == 0 ? "No files match" : $"{BatchInputs.Count} input{(BatchInputs.Count == 1 ? "" : "s")}";
    }

    public void AddBatchFiles(IEnumerable<string> files)
    {
        foreach (var f in files.Where(File.Exists))
            if (BatchInputs.All(b => b.Path != f)) BatchInputs.Add(new BatchInput { Path = f, Topology = TopologyFor(f) });
        ShortenBatchPaths();
        BatchState = $"{BatchInputs.Count} inputs";
    }

    /// <summary>The file that gives a trajectory its atoms, found beside it: a LAMMPS dump or DCD its .data; a GROMACS
    /// .xtc or .trr its .top (types, charges, bonds) or .gro (the same stem first, then the folder's only one, then
    /// topol.top / conf.gro); a .gro its own .top; an AMBER restart or trajectory its .prmtop. "" when there is none or the choice is not clear.</summary>
    internal static string TopologyFor(string f)
    {
        var ext = Path.GetExtension(f).ToLowerInvariant();
        var dir = Path.GetDirectoryName(Path.GetFullPath(f)) ?? ".";
        string Same(string e) { var p = Path.ChangeExtension(f, e); return File.Exists(p) ? p : ""; }
        string Only(string pattern) { try { var m = Directory.GetFiles(dir, pattern); return m.Length == 1 ? m[0] : ""; } catch { return ""; } }
        string Named(params string[] names) => names.Select(n => Path.Combine(dir, n)).FirstOrDefault(File.Exists) ?? "";
        string First(params string[] c) => c.FirstOrDefault(x => x.Length > 0) ?? "";
        return ext switch
        {
            ".lammpstrj" or ".dump" or ".dcd" => First(Same(".data"), Only("*.data"), ext == ".dcd" ? First(Same(".pdb"), Same(".prmtop"), Same(".parm7")) : ""),
            ".inpcrd" or ".rst7" or ".restrt" or ".rst" or ".ncrst" or ".nc" or ".mdcrd" => First(Same(".prmtop"), Same(".parm7"), Only("*.prmtop"), Only("*.parm7")),
            ".xtc" or ".trr" => First(Same(".top"), Same(".gro"), Named("topol.top"), Only("*.top"), Only("*.gro"), Named("conf.gro", "confout.gro")),
            ".gro" => Same(".top"),
            _ => "",
        };
    }

    internal static IEnumerable<string> Glob(string pattern)
    {
        if (string.IsNullOrWhiteSpace(pattern)) yield break;
        pattern = pattern.Trim();
        if (pattern.StartsWith('~')) pattern = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), pattern.TrimStart('~', '/'));
        if (!pattern.Contains('*') && !pattern.Contains('?')) { if (File.Exists(pattern)) yield return pattern; yield break; }
        var full = Path.GetFullPath(pattern);
        var root = Path.GetPathRoot(full) ?? "/";
        var parts = full[root.Length..].Split(Path.DirectorySeparatorChar, StringSplitOptions.RemoveEmptyEntries);
        IEnumerable<string> current = [root];
        for (int k = 0; k < parts.Length; ++k)
        {
            var part = parts[k];
            var last = k == parts.Length - 1;
            current = current.SelectMany(dir =>
            {
                if (!Directory.Exists(dir)) return Enumerable.Empty<string>();
                if (!part.Contains('*') && !part.Contains('?')) { var p = Path.Combine(dir, part); return last ? (File.Exists(p) ? [p] : []) : (Directory.Exists(p) ? [p] : []); }
                try { return (last ? Directory.GetFiles(dir, part) : Directory.GetDirectories(dir, part)).OrderBy(x => x, StringComparer.Ordinal); }
                catch { return Enumerable.Empty<string>(); }
            }).ToList();
        }
        foreach (var f in current) yield return f;
    }

    public void StopBatch() { _batchStop = true; BatchState = "Stopping after the inputs that are running"; }

    /// <summary>Runs the pipeline on every input (on its first or last frame), a few at a time, and writes results.csv.</summary>
    public async Task RunBatch(string? outDir = null)
    {
        if (BatchRunning || BatchInputs.Count == 0) return;
        var pipeline = PipelineRows.Count == 0 ? "" : PipelineJson();
        var sha = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(pipeline))).ToLowerInvariant();
        var last = _batchFrame == 0;
        var workers = new[] { 1, 2, 4, 8 }[_batchWorkers];
        foreach (var b in BatchInputs) { b.State = "queued"; b.Time = "—"; b.Error = ""; b.Values = []; b.Attributes.Clear(); }
        BatchRunning = true;
        _batchStop = false;
        BatchErrors = "";
        var sw = Stopwatch.StartNew();
        var queue = new Queue<BatchInput>(BatchInputs);
        // results.csv next to the inputs (batch/results.csv), or where asked; each input's pipeline outputs under outputs/
        var dir = outDir ?? Path.Combine(Path.GetDirectoryName(Path.GetDirectoryName(BatchInputs[0].Path) ?? "") ?? ".", "batch");
        var writeOutputs = PipelineOutputs.Count > 0 && pipeline.Length > 0;
        var index = BatchInputs.Select((b, k) => (b, k)).ToDictionary(x => x.b, x => x.k + 1);
        var gate = new object();
        async Task Worker()
        {
            while (true)
            {
                BatchInput? b;
                lock (gate) { if (_batchStop || queue.Count == 0) return; b = queue.Dequeue(); }
                b.State = "running";   // awaits resume on the UI thread in the app, so rows update in place
                var t0 = Stopwatch.StartNew();
                try
                {
                    var attrs = await Task.Run(() =>
                    {
                        using var doc = CapsDocument.Open(b.Path, b.Topology.Length > 0 ? b.Topology : null);
                        var frames = (int)doc.Summary().Frames;
                        if (last && frames > 1) doc.SetFrame(frames - 1);
                        doc.SetPipeline(pipeline);
                        var result = doc.PipelineResult();
                        if (writeOutputs)
                        {
                            var sub = Path.GetFileName(Path.GetDirectoryName(b.Path)) is { Length: > 0 } d ? d : Path.GetFileNameWithoutExtension(b.Path);
                            doc.PipelineWriteOutputs(pipeline, Path.Combine(dir, "outputs", $"{index[b]}_{sub}"));
                        }
                        var map = new Dictionary<string, double>();
                        if (result.Length > 0 && JsonNode.Parse(result)?["attributes"] is JsonArray a)
                            foreach (var x in a) map[(string?)x?["name"] ?? ""] = (double?)x?["value"] ?? double.NaN;
                        else
                        {
                            var s = doc.Summary();
                            map["Particles"] = s.Atoms; map["Bonds"] = s.Bonds; map["Molecules"] = s.Molecules; map["Density"] = s.Density; map["Mass"] = s.TotalMass;
                        }
                        return map;
                    });
                    foreach (var (k, v) in attrs) b.Attributes[k] = v;
                    b.State = "done";
                    b.Time = FormatSeconds(t0.Elapsed.TotalSeconds);
                    RefreshBatchTable();
                }
                catch (Exception e)
                {
                    b.State = "failed";
                    b.Error = e.Message.Replace(b.Path + ": ", "");
                    b.Time = FormatSeconds(t0.Elapsed.TotalSeconds);
                    BatchErrors = string.Join("\n", BatchInputs.Where(x => x.State == "failed").Select(x => $"{x.Shown} · {x.Error}. Skipped; the row stays as failed."));
                }
            }
        }
        await Task.WhenAll(Enumerable.Range(0, Math.Min(workers, BatchInputs.Count)).Select(_ => Worker()));
        RefreshBatchTable();
        var done = BatchInputs.Count(b => b.State == "done");
        var failed = BatchInputs.Count(b => b.State == "failed");
        try
        {
            Directory.CreateDirectory(dir);
            var csv = new StringBuilder();
            csv.AppendLine(string.Join(",", new[] { "input", "state", "seconds", "pipeline_sha256" }.Concat(_batchKeys)));
            foreach (var b in BatchInputs)
                csv.AppendLine(string.Join(",", new[] { Quote(b.Path), b.State, b.Time.Replace(" s", "").Replace("—", ""), sha }
                    .Concat(_batchKeys.Select(k => b.Attributes.TryGetValue(k, out var v) && double.IsFinite(v) ? v.ToString("R", CultureInfo.InvariantCulture) : ""))));
            File.WriteAllText(Path.Combine(dir, "results.csv"), csv.ToString());
            File.WriteAllText(Path.Combine(dir, "pipeline.json"), pipeline.Length > 0 ? pipeline : "{\"steps\":[]}");
            BatchOut = Path.Combine(dir, "results.csv");
        }
        catch (Exception e) { BatchOut = "results.csv not written: " + e.Message; }
        BatchState = $"{done} of {BatchInputs.Count} done · {failed} failed · {FormatSeconds(sw.Elapsed.TotalSeconds)}" + (_batchStop ? " · stopped" : "");
        Status = "Batch: " + BatchState;
        BatchRunning = false;
    }

    private static string Quote(string s) => s.Contains(',') || s.Contains('"') ? "\"" + s.Replace("\"", "\"\"") + "\"" : s;
    private static string FormatSeconds(double s) => s < 60 ? s.ToString("F1", CultureInfo.InvariantCulture) + " s" : (s / 60).ToString("F1", CultureInfo.InvariantCulture) + " min";

    /// <summary>Columns: the attributes the steps add, then density; values per input; mean and s.d. over the finished ones.</summary>
    private void RefreshBatchTable()
    {
        var baseKeys = new[] { "SourceFrame", "Timestep", "Particles", "Bonds", "Molecules", "CellVolume", "Mass", "TotalCharge", "Selected" };
        var keys = BatchInputs.SelectMany(b => b.Attributes.Keys).Distinct().Where(k => !baseKeys.Contains(k) && k != "Density").ToList();
        keys.Add("Density");
        keys.AddRange(new[] { "Particles", "Molecules" });
        _batchKeys = keys;
        if (!keys.SequenceEqual(BatchColumns))
        {
            BatchColumns.Clear();
            foreach (var k in keys) BatchColumns.Add(k);
            _batchMetric = Math.Clamp(_batchMetric, 0, keys.Count - 1);
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { Raise(nameof(BatchMetric)); Raise(nameof(BatchMetricName)); });
        }
        var inv = CultureInfo.InvariantCulture;
        foreach (var b in BatchInputs)
            b.Values = keys.Select(k => b.Attributes.TryGetValue(k, out var v) && double.IsFinite(v) ? (Math.Abs(v) >= 1000 || Math.Abs(v - Math.Round(v)) < 1e-9 ? v.ToString("0.##", inv) : v.ToString("G4", inv)) : "—").ToArray();
        BatchSummary.Clear();
        var done = BatchInputs.Where(b => b.State == "done").ToList();
        foreach (var k in keys.Take(8))
        {
            var xs = done.Select(b => b.Attributes.TryGetValue(k, out var v) ? v : double.NaN).Where(double.IsFinite).ToList();
            if (xs.Count == 0) continue;
            var mean = xs.Average();
            var sd = xs.Count > 1 ? Math.Sqrt(xs.Sum(x => (x - mean) * (x - mean)) / (xs.Count - 1)) : 0;
            BatchSummary.Add(new Row(k, $"{mean.ToString("G5", inv)} ± {sd.ToString("G3", inv)}  (n = {xs.Count})"));
        }
        BatchPlotChanged?.Invoke();
    }

    /// <summary>The chosen metric per input (input number, value) for the finished inputs.</summary>
    public (double X, double Y)[] BatchPlotPoints()
    {
        if (_batchMetric < 0 || _batchMetric >= _batchKeys.Count) return [];
        var k = _batchKeys[_batchMetric];
        return BatchInputs.Select((b, i) => (X: (double)(i + 1), Y: b.State == "done" && b.Attributes.TryGetValue(k, out var v) ? v : double.NaN))
                          .Where(p => double.IsFinite(p.Y)).ToArray();
    }
    /// <summary>The plotted metric by name (the picker binds to it, so it survives the columns being rebuilt).</summary>
    public string? BatchMetricName
    {
        get => _batchMetric >= 0 && _batchMetric < _batchKeys.Count ? _batchKeys[_batchMetric] : null;
        set { var i = value == null ? -1 : _batchKeys.IndexOf(value); if (i >= 0 && i != _batchMetric) { _batchMetric = i; Raise(); BatchPlotChanged?.Invoke(); } }
    }

    /// <summary>Inputs are shown relative to the folder they share (seed21/PS_melt.data).</summary>
    private void ShortenBatchPaths()
    {
        if (BatchInputs.Count < 2) return;
        var common = System.IO.Path.GetDirectoryName(BatchInputs[0].Path) ?? "";
        foreach (var b in BatchInputs)
            while (common.Length > 0 && !b.Path.StartsWith(common + System.IO.Path.DirectorySeparatorChar, StringComparison.Ordinal))
                common = System.IO.Path.GetDirectoryName(common) ?? "";
        foreach (var b in BatchInputs) b.Shown = common.Length > 0 ? b.Path[(common.Length + 1)..] : b.Path;
    }
}
