using System.Diagnostics;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Large systems (design/boards/MillionAtoms): level of detail by distance from the focus — full atoms and bonds
/// near, spheres without bonds in the middle, points far — a performance readout (frame time, atoms per tier), the memory
/// the document holds, and an adaptive mode that pulls the near radius in to keep the view under 16.7 ms.</summary>
public sealed partial class MainViewModel
{
    public const int LodAutoAtoms = 200_000;
    private bool _lodOpen, _lodOn, _lodAoNear = true, _lodAdaptive = true, _perfHud;
    private decimal _lodNear = 40, _lodFar = 80;
    private double _lodNearNow = 40, _frameMs;
    private string _lodTiers = "", _lodMemory = "", _lodBench = "";

    public bool LodOpen { get => _lodOpen; set { if (Set(ref _lodOpen, value)) { if (value) { AppearanceOpen = false; SelectionOpen = false; InteractionsOpen = false; HistoryOpen = false; LensOpen = false; RefreshLod(); } Raise(nameof(ShowLodPanel)); Raise(nameof(ShowStudioTabs)); } } }
    public bool ShowLodPanel => IsStudio && _lodOpen && _doc != null;
    public bool LodOn { get => _lodOn; set { if (Set(ref _lodOn, value)) { _lodNearNow = (double)_lodNear; RenderRequested?.Invoke(); } } }
    public decimal LodNear { get => _lodNear; set { if (Set(ref _lodNear, Math.Clamp(value, 5, 10000))) { if (_lodFar < _lodNear) LodFar = _lodNear; _lodNearNow = (double)_lodNear; RenderRequested?.Invoke(); } } }
    public decimal LodFar { get => _lodFar; set { if (Set(ref _lodFar, Math.Clamp(value, _lodNear, 20000))) RenderRequested?.Invoke(); } }
    public bool LodAoNear { get => _lodAoNear; set { if (Set(ref _lodAoNear, value)) RenderRequested?.Invoke(); } }
    public bool LodAdaptive { get => _lodAdaptive; set => Set(ref _lodAdaptive, value); }
    public bool PerfHud { get => _perfHud; set => Set(ref _perfHud, value); }
    public string LodTiers { get => _lodTiers; private set => Set(ref _lodTiers, value); }
    public string LodMemory { get => _lodMemory; private set => Set(ref _lodMemory, value); }
    public string LodBench { get => _lodBench; private set => Set(ref _lodBench, value); }
    public string FrameText => _frameMs <= 0 ? "—" : $"{_frameMs:0} ms" + (_frameMs <= 16.7 ? "" : " · above target");
    public string PerfAtoms => _doc == null ? "—" : _doc.Summary().Atoms.ToString("N0", CultureInfo.InvariantCulture);
    public string LodHudText => !_lodOn ? "off" : $"points > {_lodFar:0} Å";

    /// <summary>Level of detail in the view options (the near radius as adapted).</summary>
    private void ApplyLod(ref CapsRenderOpts o)
    {
        if (!_lodOn) return;
        o.LodNear = _lodNearNow;
        o.LodFar = Math.Max(_lodNearNow, (double)_lodFar);
        if (!_lodAoNear) o.AmbientOcclusion = 0;
    }

    /// <summary>The view's frame time: adaptive LOD pulls the near radius in above 16.7 ms and lets it back out below 8 ms.</summary>
    public void ReportFrame(double ms)
    {
        _frameMs = ms;
        Raise(nameof(FrameText));
        if (_lodOn && _lodAdaptive)
        {
            var target = (double)_lodNear;
            if (ms > 16.7 && _lodNearNow > 8) { _lodNearNow = Math.Max(8, _lodNearNow * 0.8); RenderRequested?.Invoke(); }
            else if (ms < 8 && _lodNearNow < target) { _lodNearNow = Math.Min(target, _lodNearNow * 1.25); RenderRequested?.Invoke(); }
        }
        if (_perfHud || _lodOpen) RefreshLodTiers();
    }

    private void RefreshLodTiers()
    {
        if (_doc == null) return;
        var (n, m, f, b) = _doc.RenderStats();
        LodTiers = _lodOn ? $"near {n:N0} · mid {m:N0} · points {f:N0} · {b:N0} bond halves · near radius now {_lodNearNow:0} Å" : $"all {n:N0} atoms in full · {b:N0} bond halves";
        Raise(nameof(PerfAtoms));
        Raise(nameof(LodHudText));
    }

    public void RefreshLod()
    {
        if (_doc == null) return;
        try
        {
            var j = JsonNode.Parse(_doc.Memory())!;
            var mb = (j["topology_bytes"]!.GetValue<double>() + j["frame_bytes"]!.GetValue<double>()) / 1e6;
            LodMemory = $"{j["per_atom_bytes"]!.GetValue<double>():0} B per atom ({j["atom_struct_bytes"]!.GetValue<double>():0} B atom record, 24 B per frame position) · " +
                        $"{mb:0.0} MB for {j["atoms"]!.GetValue<double>():N0} atoms × {j["frames"]!.GetValue<double>():0} frames";
        }
        catch (Exception e) { LodMemory = e.Message; }
        RefreshLodTiers();
    }

    /// <summary>Ten renders of this view at its size; mean and best frame time.</summary>
    public async Task BenchmarkView(int w, int h)
    {
        if (_doc == null) return;
        var doc = _doc;
        var opt = ViewOptions(w, h, 1);
        var cam = Camera;
        LodBench = "Rendering…";
        var times = await Task.Run(() =>
        {
            var buf = new byte[w * h * 4];
            var t = new List<double>();
            for (int k = 0; k < 10; k++)
            {
                var sw = Stopwatch.StartNew();
                doc.Render(cam, opt, buf);
                t.Add(sw.Elapsed.TotalMilliseconds);
            }
            return t;
        });
        LodBench = $"{w}×{h} px · mean {times.Average():0} ms · best {times.Min():0} ms · {1000 / times.Average():0.#} frames/s on the CPU renderer";
    }

    /// <summary>Large documents open with the level of detail on.</summary>
    private void AutoLod(CapsDocument doc)
    {
        var n = doc.Summary().Atoms;
        if (n >= LodAutoAtoms && !_lodOn) { LodOn = true; PerfHud = true; Status = $"{n:N0} atoms: level of detail on (Performance panel)"; }
    }
}
