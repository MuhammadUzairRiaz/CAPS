using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Diagnostics;
using System.Globalization;
using Avalonia.Threading;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One stage of opening a file: queued, running (with its share done), done or stopped.</summary>
public sealed class LoadStage : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    public string Title { get; init; } = "";
    private string _detail = "", _state = "queued";
    private double _fraction;
    public string Detail { get => _detail; set { _detail = value; Raise(nameof(Detail)); } }
    public string State
    {
        get => _state;
        set { _state = value; Raise(nameof(State)); Raise(nameof(Icon)); Raise(nameof(Colour)); Raise(nameof(IsRunning)); }
    }
    public double Fraction { get => _fraction; set { _fraction = value; Raise(nameof(Fraction)); Raise(nameof(Percent)); } }
    public double Percent => 100 * _fraction;
    public bool IsRunning => _state == "running";
    public string Icon => _state switch { "done" => "check", "running" => "play", "stopped" => "alert", _ => "pause" };
    public Avalonia.Media.IBrush Colour => Tokens.Brush(_state switch { "done" => "OkB", "running" => "AccB", "stopped" => "WarnB", _ => "DimB" });
}

/// <summary>Opening a large trajectory (design/boards/VisLoading): frame 0 is read and shown first, so orbit, pick and
/// measure work at once; the other frames are read in the background into the same document; then the file checks run.</summary>
public sealed partial class MainViewModel
{
    /// <summary>Dumps at least this large open progressively; smaller files open at once.</summary>
    public static long ProgressiveBytes = long.TryParse(Environment.GetEnvironmentVariable("CAPS_PROGRESSIVE_BYTES"), out var b) ? b : 4L << 20;

    public ObservableCollection<LoadStage> LoadStages { get; } = new();
    private bool _loading, _loadBackground, _loadIndexing;
    private volatile bool _loadStop;
    private string _loadTitle = "", _loadStep = "", _loadSize = "", _loadTimeLeft = "", _loadFramesText = "", _loadEndText = "";
    private double _loadFraction;
    private int _loadFramesRead;
    private CapsDocument? _loadDoc;
    private int _loadGen;   // posts from an earlier open, or arriving after reading ended, are dropped

    public bool IsLoading { get => _loading; private set { if (Set(ref _loading, value)) Raise(nameof(ShowLoadPanel)); } }
    public bool LoadInBackground { get => _loadBackground; private set { if (Set(ref _loadBackground, value)) Raise(nameof(ShowLoadPanel)); } }
    public bool ShowLoadPanel => _loading && !_loadBackground;
    /// <summary>Frame 0 is shown and the other frames are being read.</summary>
    public bool LoadIndexing { get => _loadIndexing; private set => Set(ref _loadIndexing, value); }
    public string LoadTitle { get => _loadTitle; private set => Set(ref _loadTitle, value); }
    public string LoadStep { get => _loadStep; private set => Set(ref _loadStep, value); }
    public string LoadSize { get => _loadSize; private set => Set(ref _loadSize, value); }
    public string LoadTimeLeft { get => _loadTimeLeft; private set => Set(ref _loadTimeLeft, value); }
    public double LoadFraction { get => _loadFraction; private set { if (Set(ref _loadFraction, value)) Raise(nameof(LoadPercent)); } }
    public double LoadPercent => 100 * _loadFraction;
    public string LoadFramesText { get => _loadFramesText; private set => Set(ref _loadFramesText, value); }
    public string LoadEndText { get => _loadEndText; private set => Set(ref _loadEndText, value); }
    public int LoadFramesRead { get => _loadFramesRead; private set => Set(ref _loadFramesRead, value); }

    private static readonly string[] StageTitles = ["Detect format", "Read frame 0", "Join topology", "Read the other frames", "Check the file"];

    /// <summary>True when the file should open progressively: a LAMMPS dump of at least ProgressiveBytes.</summary>
    public static bool OpensProgressively(string path)
    {
        if (!(path.EndsWith(".lammpstrj", StringComparison.OrdinalIgnoreCase) || path.EndsWith(".dump", StringComparison.OrdinalIgnoreCase))) return false;
        try { return new FileInfo(path).Length >= ProgressiveBytes; } catch { return false; }
    }

    private void SetStage(int k, string state, string? detail = null, double fraction = 0)
    {
        if (k < 0 || k >= LoadStages.Count) return;
        var s = LoadStages[k];
        s.State = state;
        if (detail != null) s.Detail = detail;
        s.Fraction = fraction;
        if (state == "running") LoadStep = $"step {k + 1} of {LoadStages.Count}";
    }

    public async Task OpenProgressive(string path, string? topology)
    {
        if (Busy) { Status = "Wait for the run to finish (or cancel it) before opening another structure"; return; }
        if (_loading) { _loadStop = true; while (_loading) await Task.Delay(20); }
        var inv = CultureInfo.InvariantCulture;
        var name = System.IO.Path.GetFileName(path);
        long size = 0;
        try { size = new FileInfo(path).Length; } catch { }
        var gen = ++_loadGen;
        LoadStages.Clear();
        foreach (var t in StageTitles) LoadStages.Add(new LoadStage { Title = t });
        _loadStop = false;
        LoadInBackground = false;
        LoadTitle = "Opening " + name;
        LoadSize = size >= 1 << 30 ? string.Format(inv, "{0:F2} GB", size / 1073741824.0) : string.Format(inv, "{0:F1} MB", size / 1048576.0);
        LoadTimeLeft = "–";
        LoadFraction = 0;
        LoadFramesRead = 0;
        LoadFramesText = "";
        LoadEndText = "";
        IsLoading = true;
        SetStage(0, "running");
        Status = LoadTitle;
        var delay = int.TryParse(Environment.GetEnvironmentVariable("CAPS_LOAD_DELAY_MS"), out var dms) ? dms : 0;   // screenshots of the loading state

        // 1. frame 0 (and the topology), shown at once
        CapsDocument doc0;
        var details = new string?[3];
        try
        {
            doc0 = await Task.Run(() => CapsDocument.OpenStaged(path, topology, 1, (st, f, d) =>
            {
                if (st < 3) details[st] = d;
                if (st < 3) Dispatcher.UIThread.Post(() => { if (gen == _loadGen && !_loadIndexing && _loading) { SetStage(st, "done", d, 1); SetStage(st + 1, "running"); } });
                return !_loadStop;
            }));
        }
        catch (Exception e)
        {
            IsLoading = false;
            Status = $"Could not open {name}: {e.Message}";
            return;
        }
        for (int k = 0; k < 3; ++k) SetStage(k, "done", details[k], 1);
        Show(doc0, name);
        if (_doc != doc0) { IsLoading = false; return; }
        _loadDoc = doc0;
        LoadIndexing = true;
        SetStage(3, "running", "reading · frame 0 is ready");
        var s0 = doc0.Summary();
        Notify(new Notice
        {
            Key = "load.stream", Severity = "info", Icon = "play", Title = "Opening a large trajectory",
            Body = $"Streaming {name} · {s0.Atoms:N0} atoms. Frames are read from disk in the background; the first frame is ready now.",
            Primary = "Cancel", OnPrimary = CancelLoad,
        });
        Status = $"Opened frame 0 of {name} · reading the other frames";

        // 2. the other frames, read in full in the background and moved into the open document
        var sw = Stopwatch.StartNew();
        int per = 0;
        CapsDocument? full = null;
        string? error = null;
        try
        {
            full = await Task.Run(() => CapsDocument.OpenStaged(path, topology, 0, (st, f, d) =>
            {
                if (st != 3) return !_loadStop;
                if (delay > 0) Thread.Sleep(delay);
                var frames = int.TryParse(d.Split(' ')[0], out var n) ? n : 0;
                if (++per % 4 == 0 || f >= 1)
                    Dispatcher.UIThread.Post(() => { if (gen == _loadGen && _loadIndexing) LoadProgress(f, frames, sw.Elapsed.TotalSeconds); });
                return !_loadStop;
            }));
        }
        catch (Exception e) { error = e.Message; }
        var stopped = _loadStop;
        LoadIndexing = false;
        Dismiss("load.stream");
        if (full != null && _doc == doc0 && _loadDoc == doc0)
        {
            var n = doc0.AdoptFrames(full);
            if (n > 0)
            {
                Frames = n;
                Raise(nameof(FrameMax));
                Raise(nameof(FrameLabel));
                RefreshSummary();
                Notes.Clear();
                foreach (var x in doc0.Notes()) Notes.Add(x);
            }
            SetStage(3, stopped ? "stopped" : "done", string.Format(inv, stopped ? "stopped · {0:N0} {1} kept" : "{0:N0} {1}", n, n == 1 ? "frame" : "frames"), 1);
        }
        else SetStage(3, "stopped", error ?? "the document was closed");
        full?.Dispose();
        _loadDoc = null;
        if (_doc != doc0) { IsLoading = false; return; }

        // 3. file checks on the whole trajectory
        SetStage(4, "running");
        LoadFileChecks();
        var look = FileChecks.Count(c => c.NeedsLook);
        SetStage(4, "done", look > 0 ? $"{look} need a look" : "nothing needs a look", 1);
        Remember(path, topology);
        IsLoading = false;
        LoadInBackground = false;
        Status = $"Opened {name} · {Frames.ToString("N0", inv)} frames" + (stopped ? " (reading stopped)" : "") + (look > 0 ? $" · {look} file check{(look == 1 ? "" : "s")} need a look" : "");
        RenderRequested?.Invoke();
    }

    private void LoadProgress(double f, int frames, double seconds)
    {
        var inv = CultureInfo.InvariantCulture;
        LoadFraction = f;
        LoadFramesRead = frames;
        var total = f > 0.02 ? (int)Math.Round(frames / f) : 0;
        LoadFramesText = total > 0 ? string.Format(inv, "{0:N0} of ~{1:N0} frames", frames, total) : string.Format(inv, "{0:N0} frames", frames);
        LoadEndText = total > 0 ? string.Format(inv, "~{0:N0}", total - 1) : "";
        SetStage(3, "running", LoadFramesText, f);
        var left = f > 0.02 ? seconds / f * (1 - f) : double.NaN;
        LoadTimeLeft = double.IsNaN(left) ? "–" : left < 60 ? string.Format(inv, "{0:F0} s", Math.Max(1, left)) : string.Format(inv, "{0:F1} min", left / 60);
        if (LoadInBackground) Status = string.Format(inv, "Reading {0} · {1:F0} % · {2}", _loadTitle.Replace("Opening ", ""), 100 * f, LoadFramesText);
    }

    /// <summary>Stops reading: the frames read so far are kept.</summary>
    public void CancelLoad()
    {
        if (!_loading) return;
        _loadStop = true;
        Status = "Stopping · the frames read so far are kept";
    }

    /// <summary>Hides the panel; reading goes on, with its progress in the status bar.</summary>
    public void LoadToBackground()
    {
        if (!_loading) return;
        LoadInBackground = true;
        Status = $"Reading {_loadTitle.Replace("Opening ", "")} in the background";
    }
}
