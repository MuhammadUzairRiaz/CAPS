using System;
using System.Linq;
using System.Text.Json.Nodes;
using Avalonia.Threading;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Motion (design/boards/Motion): camera moves fly over 450 ms with the standard easing, or cut when motion is
/// reduced (Settings, or the operating system's setting).</summary>
public partial class MainViewModel
{
    private DispatcherTimer? _fly;
    private CapsCamera _flyFrom, _flyTo;
    private DateTime _flyStart;

    public string SetReduceMotion
    {
        get => _settings.ReduceMotion;
        set
        {
            if (_settings.ReduceMotion == value) return;
            _settings.ReduceMotion = value;
            Motion.Mode = value;
            Raise(); Raise(nameof(MotionText));
            Changed("Reduce motion");
        }
    }
    public string MotionText => Motion.Reduced
        ? $"Camera moves are cuts; progress bars stay; nothing loops ({(Motion.Mode == "on" ? "reduced here" : Motion.SystemText)})."
        : $"Camera moves fly for 450 ms with the standard easing ({(Motion.Mode == "off" ? "full motion here" : Motion.SystemText)}).";

    /// <summary>Only the camera changed (a flight): the view turns what it has, nothing is rebuilt.</summary>
    public event Action? ViewRequested;

    /// <summary>Moves the camera to target: a 450 ms flight, or a cut when motion is reduced. Returns at once.</summary>
    public void FlyTo(CapsCamera target)
    {
        _fly?.Stop();
        if (Motion.Reduced || _doc == null)
        {
            Camera = target;
            RenderRequested?.Invoke();
            return;
        }
        _flyFrom = Camera;
        _flyTo = target;
        // turn the short way round
        _flyTo.Yaw = _flyFrom.Yaw + Math.IEEERemainder(target.Yaw - _flyFrom.Yaw, 2 * Math.PI);
        _flyStart = DateTime.UtcNow;
        _fly ??= new DispatcherTimer(TimeSpan.FromMilliseconds(16), DispatcherPriority.Render, (_, _) => FlyStep());
        _fly.Start();
    }

    public void StopFly() => _fly?.Stop();

    /// <summary>One frame of the flight (also called directly by tests with a time).</summary>
    public bool FlyStep(double? at = null)
    {
        var t = at ?? (DateTime.UtcNow - _flyStart).TotalMilliseconds / Motion.Camera.TotalMilliseconds;
        var done = t >= 1;
        var e = Motion.Standard(Math.Min(1, t));
        static double Lerp(double a, double b, double f) => a + (b - a) * f;
        var c = _flyTo;
        c.Yaw = Lerp(_flyFrom.Yaw, _flyTo.Yaw, e);
        c.Pitch = Lerp(_flyFrom.Pitch, _flyTo.Pitch, e);
        c.Zoom = Math.Exp(Lerp(Math.Log(Math.Max(1e-3, _flyFrom.Zoom)), Math.Log(Math.Max(1e-3, _flyTo.Zoom)), e));
        c.PanX = Lerp(_flyFrom.PanX, _flyTo.PanX, e);
        c.PanY = Lerp(_flyFrom.PanY, _flyTo.PanY, e);
        Camera = c;
        if (done) _fly?.Stop();
        if (ViewRequested != null) ViewRequested.Invoke(); else RenderRequested?.Invoke();   // only the camera moved
        return done;
    }

    /// <summary>F: the selection (named selection, else picked atoms) centred and filling 60 % of the view; nothing
    /// selected frames the whole structure.</summary>
    public void FrameSelection()
    {
        if (_doc == null) return;
        int[] idx = [];
        try
        {
            if (JsonNode.Parse(_doc.SelectionJson())?["indices"] is JsonArray a && a.Count > 0) idx = a.Select(x => (int)x!).ToArray();
        }
        catch { /* no selection support */ }
        if (idx.Length == 0) idx = _selection.ToArray();
        var target = idx.Length == 0 ? new CapsCamera { Yaw = Camera.Yaw, Pitch = Camera.Pitch, Zoom = 1, Perspective = Camera.Perspective } : _doc.Focus(Camera, idx);
        FlyTo(target);
        Status = idx.Length == 0 ? "Framed the whole structure" : $"Framed {idx.Length} selected atom{(idx.Length == 1 ? "" : "s")}";
    }
}
