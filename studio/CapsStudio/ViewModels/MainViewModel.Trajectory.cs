using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A plot of the trajectory player: the series it shows, its points and the value at the current frame.</summary>
public sealed class LinkedPlot : ObservableObject
{
    public LinkedPlot(string colour) { Colour = Avalonia.Media.Brush.Parse(colour); }
    public Avalonia.Media.IBrush Colour { get; }
    private string? _series;
    public string? Series { get => _series; set { if (value != null && Set(ref _series, value)) Changed?.Invoke(); } }
    public (double X, double Y)[] Points { get; set; } = [];
    private string _value = "—";
    public string Value { get => _value; set => Set(ref _value, value); }
    public event Action? Changed;
}

/// <summary>Trajectory player (design/boards/Trajectory): playback with speed, loop and smoothing, run starts marked on
/// the timeline, one chain's end-to-end vector drawn in the view, and plots linked to the frame (density, volume, the
/// chain's Rg and Ree, and a LAMMPS log's thermo columns matched by timestep).</summary>
public sealed partial class MainViewModel
{
    public bool IsTrajectory => _module == 32;
    public ObservableCollection<string> TrajSeriesNames { get; } = new();
    public LinkedPlot[] TrajPlots { get; } = [new("#F5A524"), new("#6CC4D8"), new("#E07A5F")];
    public static readonly string[] TrajSpeeds = ["0.5×", "1×", "2×", "4×"];
    public event Action? TrajectoryChanged;

    private JsonNode? _traj;
    private int _trajTicket;
    private string _trajLog = "", _trajNote = "", _trajError = "";
    private decimal _trajDt = 1;
    private int _trajMolecule, _trajSpeed = 1;
    private bool _trajSmooth;

    public void OpenTrajectory()
    {
        if (_doc == null) return;
        foreach (var p in TrajPlots) p.Changed -= TrajPlotChanged;
        foreach (var p in TrajPlots) p.Changed += TrajPlotChanged;
        // a LAMMPS log beside the trajectory
        if (_trajLog.Length == 0 && Path.GetDirectoryName(_doc.Path) is { Length: > 0 } dir && Directory.Exists(dir))
        {
            var log = Directory.EnumerateFiles(dir).FirstOrDefault(f => Path.GetFileName(f).StartsWith("log.", StringComparison.OrdinalIgnoreCase) || f.EndsWith(".log", StringComparison.OrdinalIgnoreCase));
            if (log != null) _trajLog = log;
        }
        SetModule(32);
        Raise(nameof(TrajLogText));
        TrajCompute();
    }

    public string TrajLogPath { get => _trajLog; set { if (Set(ref _trajLog, value ?? "")) { Raise(nameof(TrajLogText)); TrajCompute(); } } }
    public string TrajLogText => _trajLog.Length == 0 ? "no LAMMPS log" : "log: " + Path.GetFileName(_trajLog);
    public decimal TrajDt { get => _trajDt; set { if (Set(ref _trajDt, Math.Clamp(value, 0.01m, 100m))) TrajCompute(); } }
    /// <summary>The chain for Rg, Ree and the end-to-end line (0: the largest).</summary>
    public int TrajMolecule { get => _trajMolecule; set { if (Set(ref _trajMolecule, Math.Max(0, value))) TrajCompute(); } }
    public int TrajSpeed { get => _trajSpeed; set { if (Set(ref _trajSpeed, Math.Clamp(value, 0, 3))) PlaybackFps = _trajSpeed switch { 0 => 4, 1 => 8, 2 => 16, _ => 32 }; } }
    public bool TrajSmooth
    {
        get => _trajSmooth;
        set
        {
            if (!Set(ref _trajSmooth, value) || _doc == null) return;
            _doc.SetSmoothing(value ? 5 : 1);
            RenderRequested?.Invoke();
            TrajectoryChanged?.Invoke();
        }
    }
    public string TrajError { get => _trajError; private set { if (Set(ref _trajError, value)) Raise(nameof(TrajHasError)); } }
    public bool TrajHasError => _trajError.Length > 0;
    public string TrajNote { get => _trajNote; private set => Set(ref _trajNote, value); }
    public int[] TrajEnds { get; private set; } = [-1, -1];
    public List<int> TrajRunFrames { get; private set; } = new();
    public int TrajChain { get; private set; }

    /// <summary>The per-frame series (in the background; frames are read from memory).</summary>
    public void TrajCompute()
    {
        if (_doc == null || _module != 32) return;
        var doc = _doc;
        var ticket = ++_trajTicket;
        var opts = new JsonObject { ["molecule"] = _trajMolecule, ["dt_fs"] = (double)_trajDt, ["log"] = _trajLog, ["stride"] = 1 }.ToJsonString();
        TrajNote = "Computing the series…";
        Task.Run(() => doc.TrajectorySeries(opts)).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            if (ticket != _trajTicket) return;
            if (t.Exception != null) { TrajError = t.Exception.InnerException?.Message ?? "cannot compute the series"; return; }
            var j = JsonNode.Parse(t.Result)!;
            if (j["ok"]?.GetValue<bool>() != true) { TrajError = j["error"]?.GetValue<string>() ?? "cannot compute the series"; TrajNote = ""; return; }
            TrajError = "";
            _traj = j;
            TrajChain = (int)j["molecule"]!.GetValue<double>();
            TrajEnds = j["ends"]!.AsArray().Select(x => (int)x!.GetValue<double>()).ToArray();
            TrajRunFrames = j["run_frames"]!.AsArray().Select(x => (int)x!.GetValue<double>()).ToList();
            var names = j["columns"]!.AsArray().Select(x => x!.GetValue<string>()).Skip(3).ToList();   // Frame, Timestep, Time are the x axes
            TrajSeriesNames.Clear();
            foreach (var n in names) TrajSeriesNames.Add(n);
            string? Pick(params string[] prefer) => prefer.FirstOrDefault(p => names.Contains(p));
            TrajPlots[0].Series ??= Pick("Density (g/cm³)", "Volume (Å³)", "Rg (Å)");
            TrajPlots[1].Series ??= Pick("Temp", "Press", "PotEng", "Volume (Å³)", "Ree (Å)");
            TrajPlots[2].Series ??= Pick("Rg (Å)");
            foreach (var p in TrajPlots) if (p.Series != null && !names.Contains(p.Series)) p.Series = names.FirstOrDefault();
            var logRows = (int)(j["log_rows"]?.GetValue<double>() ?? 0);
            TrajNote = $"Series of {_frames:N0} frames · chain {TrajChain}" + (logRows > 0 ? $" · {logRows:N0} thermo rows from {Path.GetFileName(_trajLog)}" : "");
            Raise(nameof(TrajChainText));
            Raise(nameof(TrajHudTitle));
            Raise(nameof(TrajStatus));
            TrajPlotChanged();
        }));
    }

    public string TrajChainText => $"chain {TrajChain}";
    public string PlayIcon => IsPlaying ? "pause" : "play";

    private void TrajPlotChanged()
    {
        if (_traj == null) return;
        var cols = _traj["columns"]!.AsArray().Select(x => x!.GetValue<string>()).ToList();
        var rows = _traj["rows"]!.AsArray();
        foreach (var p in TrajPlots)
        {
            var k = p.Series == null ? -1 : cols.IndexOf(p.Series);
            p.Points = k < 0 ? [] : rows.Select(r => ((double?)r?[2] ?? 0, (double?)r?[k] ?? double.NaN)).Where(q => double.IsFinite(q.Item2)).ToArray();
        }
        TrajUpdateValues();
        TrajectoryChanged?.Invoke();
    }

    /// <summary>The values at the current frame (the plots' chips) and the time.</summary>
    public void TrajUpdateValues()
    {
        Raise(nameof(TrajFrameText));
        Raise(nameof(TrajTimeText));
        Raise(nameof(TrajCursor));
        Raise(nameof(TrajReeText));
        if (_traj == null) return;
        var cols = _traj["columns"]!.AsArray().Select(x => x!.GetValue<string>()).ToList();
        var rows = _traj["rows"]!.AsArray();
        var row = _frame < rows.Count ? rows[_frame] : null;
        foreach (var p in TrajPlots)
        {
            var k = p.Series == null ? -1 : cols.IndexOf(p.Series);
            var v = k < 0 || row == null ? null : (double?)row[k];
            p.Value = v is double d ? (Math.Abs(d) >= 1000 ? d.ToString("0.0", CultureInfo.InvariantCulture) : Math.Abs(d) >= 10 ? d.ToString("0.00", CultureInfo.InvariantCulture) : d.ToString("0.000", CultureInfo.InvariantCulture)) : "—";
        }
    }

    private double? RowValue(string column)
    {
        if (_traj == null) return null;
        var k = _traj["columns"]!.AsArray().Select(x => x!.GetValue<string>()).ToList().IndexOf(column);
        var rows = _traj["rows"]!.AsArray();
        return k < 0 || _frame >= rows.Count ? null : (double?)rows[_frame]?[k];
    }

    public string TrajFrameText => $"frame {_frame} / {Math.Max(0, _frames - 1)}";
    public string TrajTimeText => RowValue("Time (ps)") is double t ? (t >= 1000 ? $"t = {t / 1000:0.000} ns" : $"t = {t:0.00} ps") : "";
    public double? TrajCursor => RowValue("Time (ps)");
    public string TrajReeText => RowValue("Ree (Å)") is double r ? $"R_ee chain {TrajChain} · {r:0.0} Å" : "";
    public string TrajHudTitle => $"{(_doc != null ? Path.GetFileName(_doc.Path) : Title)} · {_frames:N0} frames";
    public string TrajStatus
    {
        get
        {
            var s = _doc?.Summary();
            var dtps = _traj?["rows"] is JsonArray rows && rows.Count > 1 ? ((double?)rows[1]?[2] ?? 0) - ((double?)rows[0]?[2] ?? 0) : 0;
            return $"{s?.Atoms ?? 0:N0} atoms · {_frames:N0} frames" + (dtps > 0 ? $" · {dtps:0.###} ps/frame" : "");
        }
    }
}
