using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Runtime.InteropServices;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One benchmark table: what it measures, its last result.</summary>
public sealed class BenchItem : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));

    public string Id { get; init; } = "";
    public string Title { get; init; } = "";
    public string Scope { get; init; } = "";
    public bool Runnable { get; init; } = true;
    private JsonObject _json = new();
    /// <summary>The table as caps_bench_run returns it.</summary>
    public JsonObject Json
    {
        get => _json;
        set
        {
            _json = value;
            foreach (var n in new[] { nameof(Status), nameof(RowsText), nameof(IsPass), nameof(IsFail), nameof(IsRunning), nameof(Note), nameof(Columns), nameof(Rows), nameof(SecondsText),
                                      nameof(HasFailingRow), nameof(FailingTitle), nameof(FailingFacts), nameof(FailingFile), nameof(CanOpenFailing) }) Raise(n);
        }
    }
    private bool _running;
    public bool Running { get => _running; set { _running = value; Raise(nameof(Running)); Raise(nameof(Status)); Raise(nameof(IsRunning)); } }
    public string Status => _running ? "running" : (string?)_json["status"] ?? "not run";
    public bool IsPass => Status == "pass";
    public bool IsFail => Status == "fail";
    public bool IsRunning => _running;
    public string Note => (string?)_json["note"] ?? "";
    public string SecondsText => _json["seconds"] is JsonNode s && (double)s > 0 ? $"{(double)s:F1} s" : "";
    public string[] Columns => _json["columns"] is JsonArray a ? a.Select(x => (string?)x ?? "").ToArray() : [];
    public List<(string[] Cells, string Status)> Rows =>
        _json["rows"] is JsonArray a ? a.Select(r => (((JsonArray)r!["cells"]!).Select(c => (string?)c ?? "").ToArray(), (string?)r["status"] ?? "")).ToList() : [];
    // the failing row (design/boards/Bench "Table 7 · failing row"): its cells by column and the run it came from
    private JsonObject? FirstFail => _json["rows"] is JsonArray a ? a.OfType<JsonObject>().FirstOrDefault(r => (string?)r["status"] == "fail") : null;
    public bool HasFailingRow => !_running && FirstFail != null;
    public string FailingTitle
    {
        get
        {
            var n = _json["rows"] is JsonArray a ? a.OfType<JsonObject>().Count(r => (string?)r["status"] == "fail") : 0;
            return $"{Id} · failing row" + (n > 1 ? $" (1 of {n})" : "");
        }
    }
    public List<JobFact> FailingFacts
    {
        get
        {
            if (FirstFail is not { } r || r["cells"] is not JsonArray cells) return [];
            var cols = Columns;
            return cells.Select((c, i) => new JobFact(i < cols.Length ? cols[i] : $"column {i + 1}", (string?)c ?? ""))
                        .Where(f => f.Key != "Status").ToList();
        }
    }
    public string FailingFile => (string?)FirstFail?["file"] ?? "";
    public bool CanOpenFailing => FailingFile.Length > 0 && File.Exists(FailingFile);
    public string RowsText
    {
        get
        {
            var rows = Rows;
            if (!Runnable) return "—";
            if (rows.Count == 0) return "0";
            var judged = rows.Count(r => r.Status is "pass" or "fail");
            return judged == 0 ? rows.Count.ToString() : $"{rows.Count(r => r.Status == "pass")} / {judged}";
        }
    }
}

/// <summary>Bench (design/boards/Bench): the validation suite that regenerates the paper tables.</summary>
public sealed partial class MainViewModel
{
    public bool IsBench => _module == 12;
    public ObservableCollection<BenchItem> BenchItems { get; } = new();
    private BenchItem? _bench;
    public BenchItem? SelectedBench { get => _bench; set { if (Set(ref _bench, value)) BenchSelected?.Invoke(); } }
    public event Action? BenchSelected;
    private bool _benchRunning, _benchCancel;
    public bool BenchRunning { get => _benchRunning; private set { if (Set(ref _benchRunning, value)) { Raise(nameof(BenchIdle)); } } }
    public bool BenchIdle => !_benchRunning;
    private string _benchProgress = "";
    public string BenchProgress { get => _benchProgress; private set => Set(ref _benchProgress, value); }
    private decimal _benchRepeats = 3;
    public decimal BenchRepeats { get => _benchRepeats; set { if (Set(ref _benchRepeats, Math.Clamp(value, 1, 10))) Raise(nameof(BenchCommand)); } }
    private bool _benchQuick;
    public bool BenchQuick { get => _benchQuick; set { if (Set(ref _benchQuick, value)) Raise(nameof(BenchCommand)); } }
    public string BenchCommand => $"caps bench --all --repeats {(int)_benchRepeats}{(_benchQuick ? " --quick" : "")}";
    public string BenchPassText => $"{BenchItems.Count(b => b.IsPass)} pass";
    public string BenchFailText => $"{BenchItems.Count(b => b.IsFail)} fail";
    public bool BenchAnyFail => BenchItems.Any(b => b.IsFail);
    public List<JobFact> BenchEnvironment { get; } = new();

    private static string BenchFile => AppSettings.Override != null ? Path.Combine(Path.GetDirectoryName(AppSettings.Override)!, "caps-bench.json") : Path.Combine(AppSettings.Folder, "bench.json");

    public void LoadBench()
    {
        if (BenchItems.Count > 0) return;
        try
        {
            var n = Native.BenchList(null, 0);
            var buf = new byte[Math.Max(1, n)];
            Native.BenchList(buf, buf.Length);
            var list = JsonNode.Parse(System.Text.Encoding.UTF8.GetString(buf, 0, Math.Max(0, n - 1))) as JsonArray ?? [];
            JsonArray? saved = null;
            try { if (File.Exists(BenchFile)) saved = JsonNode.Parse(File.ReadAllText(BenchFile)) as JsonArray; } catch { }
            foreach (var t in list.OfType<JsonObject>())
            {
                var id = (string?)t["id"] ?? "";
                var item = new BenchItem
                {
                    Id = id, Title = (string?)t["title"] ?? "", Scope = (string?)t["scope"] ?? "",
                    Runnable = ((string?)t["note"] ?? "").Length == 0,
                };
                var last = saved?.OfType<JsonObject>().FirstOrDefault(x => (string?)x["id"] == id);
                item.Json = (JsonObject)(last ?? t).DeepClone();
                BenchItems.Add(item);
            }
        }
        catch (Exception e) { Status = "Bench: " + e.Message; }
        SelectedBench = BenchItems.FirstOrDefault();
        BenchEnvironment.Clear();
        BenchEnvironment.AddRange(BenchFacts());
        RaiseBench();
    }

    private void RaiseBench() { Raise(nameof(BenchPassText)); Raise(nameof(BenchFailText)); Raise(nameof(BenchAnyFail)); }

    private static List<JobFact> BenchFacts()
    {
        string Tool(string name)
        {
            foreach (var dir in (System.Environment.GetEnvironmentVariable("PATH") ?? "").Split(Path.PathSeparator))
            {
                foreach (var ext in RuntimeInformation.IsOSPlatform(OSPlatform.Windows) ? new[] { ".exe", "" } : [""])
                {
                    var p = Path.Combine(dir, name + ext);
                    if (File.Exists(p)) return RecentFiles.Tilde(p);
                }
            }
            return "not found";
        }
        return
        [
            new("CAPS", $"{typeof(MainViewModel).Assembly.GetName().Version?.ToString(3)} · ABI {Native.AbiVersion()}"),
            new("OS", RuntimeInformation.OSDescription),
            new("CPU", $"{RuntimeInformation.ProcessArchitecture} · {System.Environment.ProcessorCount} threads"),
            new("precision", "double"),
            new(".NET", RuntimeInformation.FrameworkDescription),
            new("packmol", Tool("packmol")),
            new("LAMMPS", Tool("lmp")),
        ];
    }

    /// <summary>Runs one table, or every table the suite can run.</summary>
    public async Task RunBench(bool all)
    {
        if (_benchRunning) return;
        var samples = Paths.Samples;
        if (samples == null) { Status = "Bench needs the samples folder (ps_melt.data, water.pdb)"; return; }
        var targets = all ? BenchItems.Where(b => b.Runnable).ToList() : _bench is { Runnable: true } b1 ? [b1] : [];
        if (targets.Count == 0) { Status = "This table is not part of the built-in suite"; return; }
        BenchRunning = true;
        _benchCancel = false;
        var repeats = (int)_benchRepeats;
        var quick = _benchQuick;
        var ff = Paths.ForceFields;
        try
        {
            foreach (var t in targets)
            {
                if (_benchCancel) break;
                t.Running = true;
                SelectedBench = t;
                var json = await Task.Run(() => Native.BenchRunJson(t.Id, samples, ff, repeats, quick, (id, what, f) =>
                {
                    Avalonia.Threading.Dispatcher.UIThread.Post(() => BenchProgress = $"{id} · {what} · {f * 100:F0} %");
                    return !_benchCancel;
                }));
                t.Running = false;
                if (JsonNode.Parse(json) is JsonObject o) t.Json = o;
                RaiseBench();
                BenchSelected?.Invoke();
            }
            SaveBench();
            BenchProgress = _benchCancel ? "Cancelled" : $"Finished · {BenchPassText} · {BenchFailText}";
            Status = "Bench: " + BenchProgress;
        }
        catch (Exception e) { BenchProgress = "Could not run: " + e.Message; }
        finally
        {
            foreach (var t in targets) t.Running = false;
            BenchRunning = false;
        }
    }

    public void CancelBench() => _benchCancel = true;

    private void SaveBench()
    {
        try
        {
            Directory.CreateDirectory(Path.GetDirectoryName(BenchFile)!);
            File.WriteAllText(BenchFile, new JsonArray(BenchItems.Select(b => (JsonNode)b.Json.DeepClone()).ToArray()).ToJsonString());
        }
        catch { }
    }

    public void ExportBench(string dir)
    {
        var json = new JsonArray(BenchItems.Select(b => (JsonNode)b.Json.DeepClone()).ToArray()).ToJsonString();
        if (Native.BenchWrite(json, dir) < 0) throw new InvalidOperationException(Native.LastError());
        Status = $"Wrote results.md, results.tex and one CSV per table to {dir}";
    }
}
