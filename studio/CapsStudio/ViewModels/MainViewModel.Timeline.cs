using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Timeline (design/boards/Timeline): the pipeline at this frame, playback (every n-th frame, frames per
/// second, loop) and, once a time series is computed, an attribute's sparkline over the frames with markers where
/// another attribute changes (cluster count by default).</summary>
public sealed partial class MainViewModel
{
    private int _playFps = 8, _playEvery = 1;
    private bool _playLoop = true, _showMarkers = true;
    private string? _sparkAttr, _markerAttr;
    public event Action? PlaybackChanged;
    public event Action? TimelineChanged;
    public ObservableCollection<string> TimelineAttributes { get; } = new();

    public int PlaybackFps { get => _playFps; set { if (Set(ref _playFps, Math.Clamp(value, 1, 60))) PlaybackChanged?.Invoke(); } }
    public int PlaybackEvery { get => _playEvery; set => Set(ref _playEvery, Math.Clamp(value, 1, 10000)); }
    public bool PlaybackLoop { get => _playLoop; set => Set(ref _playLoop, value); }
    public bool ShowMarkers { get => _showMarkers; set { if (Set(ref _showMarkers, value)) TimelineChanged?.Invoke(); } }
    public string? SparkAttribute { get => _sparkAttr; set { if (value != null && Set(ref _sparkAttr, value)) TimelineChanged?.Invoke(); } }
    public string? MarkerAttribute { get => _markerAttr; set { if (value != null && Set(ref _markerAttr, value)) { Raise(nameof(MarkerNote)); TimelineChanged?.Invoke(); } } }
    public bool HasTimeline => IsVisualize && _series?["rows"] is JsonArray { Count: > 1 };
    public string MarkerNote => MarkerFrames().Count is var n && HasTimeline ? $"{n} change{(n == 1 ? "" : "s")} in {SeriesFrames().Length} frames" : "compute the time series (Data tables › Over all frames)";
    public string PipelineAtFrameText => PipelineRows.Count == 0 ? "No steps." :
        string.Join("\n", PipelineRows.Reverse().Select(r => $"{(r.Level == "error" ? "✕" : r.Enabled ? "✓" : "·")} {r.Title} — {r.Summary}"));

    /// <summary>One playback tick: forward by the stride, looping or stopping at the end.</summary>
    public bool PlayStep()
    {
        if (_doc == null || _frames < 2) return false;
        var next = _frame + _playEvery;
        if (next > _frames - 1)
        {
            if (!_playLoop) { IsPlaying = false; return false; }
            next = 0;
        }
        Frame = next;
        return true;
    }

    private int ColumnOf(string? name) => name == null || _series?["columns"] is not JsonArray c ? -1 : c.Select(x => (string?)x).ToList().IndexOf(name);

    public int[] SeriesFrames() => _series?["rows"] is JsonArray rows ? rows.Select(r => (int)((double?)r?[0] ?? 0)).ToArray() : [];

    /// <summary>The sparkline: (frame, value) of the chosen attribute.</summary>
    public (double X, double Y)[] Sparkline()
    {
        var k = ColumnOf(_sparkAttr);
        if (k < 0 || _series?["rows"] is not JsonArray rows) return [];
        return rows.Select(r => ((double?)r?[0] ?? 0, (double?)r?[k] ?? double.NaN)).Where(p => double.IsFinite(p.Item2)).ToArray();
    }

    private int _seriesWindow = 1;
    public static readonly int[] SeriesWindows = [1, 3, 5, 9, 15, 25, 51];
    /// <summary>Running mean over this many computed frames, centred (1: off), drawn over the sparkline.</summary>
    public int SeriesWindow { get => _seriesWindow; set { if (Set(ref _seriesWindow, Math.Max(1, value))) { TimelineChanged?.Invoke(); if (_series != null) LoadPipeTable(); } } }

    /// <summary>The sparkline's centred running mean over SeriesWindow points (shorter at the ends); empty when off.</summary>
    public (double X, double Y)[] SparklineMean()
    {
        var pts = Sparkline();
        if (_seriesWindow <= 1 || pts.Length < 3) return [];
        var half = _seriesWindow / 2;
        var outp = new (double X, double Y)[pts.Length];
        for (var i = 0; i < pts.Length; ++i)
        {
            int a = Math.Max(0, i - half), b = Math.Min(pts.Length - 1, i + half);
            double sum = 0;
            for (var k = a; k <= b; ++k) sum += pts[k].Y;
            outp[i] = (pts[i].X, sum / (b - a + 1));
        }
        return outp;
    }

    /// <summary>Frames where the marker attribute changes value.</summary>
    public List<int> MarkerFrames()
    {
        var k = ColumnOf(_markerAttr);
        var list = new List<int>();
        if (k < 0 || _series?["rows"] is not JsonArray rows) return list;
        double? prev = null;
        foreach (var r in rows)
        {
            var v = (double?)r?[k];
            if (prev != null && v != null && Math.Abs(v.Value - prev.Value) > 1e-9) list.Add((int)((double?)r?[0] ?? 0));
            prev = v;
        }
        return list;
    }

    /// <summary>Called when a time series arrives: offers its attributes and picks sensible defaults.</summary>
    private void RefreshTimeline()
    {
        TimelineAttributes.Clear();
        if (_series?["columns"] is JsonArray cols)
            foreach (var c in cols.Skip(1)) if ((string?)c is { } n && n != "Timestep") TimelineAttributes.Add(n);
        var skip = new[] { "Particles", "Bonds", "Molecules", "CellVolume", "Mass", "TotalCharge", "Selected" };
        if (_sparkAttr == null || !TimelineAttributes.Contains(_sparkAttr))
            _sparkAttr = TimelineAttributes.FirstOrDefault(a => !skip.Contains(a) && a != "Density" && !a.StartsWith("ClusterAnalysis")) ?? TimelineAttributes.FirstOrDefault(a => a == "Density");
        if (_markerAttr == null || !TimelineAttributes.Contains(_markerAttr))
            _markerAttr = TimelineAttributes.FirstOrDefault(a => a == "ClusterAnalysis.cluster_count") ?? TimelineAttributes.FirstOrDefault(a => a.EndsWith("count"));
        Avalonia.Threading.Dispatcher.UIThread.Post(() => { Raise(nameof(SparkAttribute)); Raise(nameof(MarkerAttribute)); });
        Raise(nameof(HasTimeline));
        Raise(nameof(MarkerNote));
        TimelineChanged?.Invoke();
    }
}
