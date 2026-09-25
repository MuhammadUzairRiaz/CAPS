using System.Collections.ObjectModel;
using System.ComponentModel;
using Avalonia.Media.Imaging;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One of the four views: its camera, picture and whether it is the active one.</summary>
public sealed class ViewportTile : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    public string Name { get; init; } = "";
    public bool Ortho { get; init; }
    public double Yaw { get; set; }
    public double Pitch { get; set; }
    private Bitmap? _image;
    private bool _active, _shown = true;
    /// <summary>Where the tile sits in the page, as fractions of its width and height.</summary>
    public double NX { get; set; }
    public double NY { get; set; }
    public double NW { get; set; } = 0.5;
    public double NH { get; set; } = 0.5;
    public Bitmap? Image { get => _image; set { _image = value; Raise(nameof(Image)); } }
    public bool Active { get => _active; set { _active = value; Raise(nameof(Active)); } }
    public bool Shown { get => _shown; set { _shown = value; Raise(nameof(Shown)); } }
    public string Label => Name + (Ortho ? " · ortho" : "");
    public int PixelW { get; set; } = 400;
    public int PixelH { get; set; } = 300;
    internal int Gen;
}

/// <summary>Studio › Viewports (design/boards/FourViewports): top, front and left orthographic views and the perspective
/// view; the ortho views share zoom and pan; a pick in any view selects in all.</summary>
public sealed partial class MainViewModel
{
    public bool IsViewports => _module == 25;
    public static readonly string[] ViewportLayouts = ["1", "2 × 1", "2 × 2", "1 + 3"];
    public ObservableCollection<ViewportTile> ViewportTiles { get; } =
    [
        new ViewportTile { Name = "Top", Ortho = true, Yaw = 0, Pitch = Math.PI / 2 },
        new ViewportTile { Name = "Front", Ortho = true, Yaw = 0, Pitch = 0 },
        new ViewportTile { Name = "Left", Ortho = true, Yaw = Math.PI / 2, Pitch = 0 },
        new ViewportTile { Name = "Perspective", Ortho = false },
    ];
    private int _vpLayout = 2, _vpActive = 3;
    private double _orthoZoom = 1, _orthoPanX, _orthoPanY;
    private bool _vpLinkZoom = true;

    public void OpenViewports()
    {
        if (_doc == null) { Status = "Open a structure first"; return; }
        SetModule(25);
        ArrangeViewports();
        RenderViewports();
    }

    public int ViewportLayout { get => _vpLayout; set { if (Set(ref _vpLayout, Math.Clamp(value, 0, 3))) { ArrangeViewports(); RenderViewports(); } } }
    public bool ViewportLinkZoom { get => _vpLinkZoom; set => Set(ref _vpLinkZoom, value); }
    public string ViewportActiveName => ViewportTiles[_vpActive].Name;
    public ViewportTile ActiveViewport => ViewportTiles[_vpActive];

    /// <summary>Places the tiles: the active one alone (1), beside a second (2 × 1), all four (2 × 2), or large with the
    /// other three stacked beside it (1 + 3).</summary>
    private void ArrangeViewports()
    {
        var active = ViewportTiles[_vpActive];
        var others = ViewportTiles.Where(t => t != active).ToList();
        foreach (var t in ViewportTiles) { t.Shown = false; t.Active = t == active; }
        void Place(ViewportTile t, double x, double y, double w, double h) { t.Shown = true; t.NX = x; t.NY = y; t.NW = w; t.NH = h; }
        switch (_vpLayout)
        {
            case 0: Place(active, 0, 0, 1, 1); break;
            case 1: Place(active, 0, 0, 0.5, 1); Place(_vpActive == 3 ? ViewportTiles[0] : ViewportTiles[3], 0.5, 0, 0.5, 1); break;
            case 3:
                Place(active, 0, 0, 0.66, 1);
                for (int k = 0; k < 3; ++k) Place(others[k], 0.66, k / 3.0, 0.34, 1 / 3.0);
                break;
            default: for (int k = 0; k < 4; ++k) Place(ViewportTiles[k], k % 2 * 0.5, k / 2 * 0.5, 0.5, 0.5); break;
        }
        Raise(nameof(ViewportActiveName));
        Raise(nameof(ActiveViewport));
        ViewportsArranged?.Invoke();
    }
    public event Action? ViewportsArranged;

    public void ActivateViewport(ViewportTile t)
    {
        var i = ViewportTiles.IndexOf(t);
        if (i < 0 || i == _vpActive) return;
        _vpActive = i;
        foreach (var x in ViewportTiles) x.Active = x == t;
        if (_vpLayout is 0 or 1 or 3) ArrangeViewports();
        Raise(nameof(ViewportActiveName));
        Raise(nameof(ActiveViewport));
    }
    public void CycleViewport() => ActivateViewport(ViewportTiles[(_vpActive + 1) % 4]);
    public void MaximiseViewport() => ViewportLayout = _vpLayout == 0 ? 2 : 0;

    private CapsCamera CameraOf(ViewportTile t)
    {
        if (!t.Ortho) { var c = Camera; c.Perspective = 1; return c; }
        return new CapsCamera { Yaw = t.Yaw, Pitch = t.Pitch, Zoom = _orthoZoom, PanX = _orthoPanX, PanY = _orthoPanY, Perspective = 0 };
    }

    /// <summary>Drag: the perspective view orbits; an ortho view pans (all ortho views when zoom is linked).</summary>
    public void DragViewport(ViewportTile t, double dx, double dy, double tileHeight)
    {
        if (!t.Ortho)
        {
            Camera.Yaw += dx * 0.008;
            Camera.Pitch = Math.Clamp(Camera.Pitch + dy * 0.008, -Math.PI / 2, Math.PI / 2);
            RenderViewport(t);
            return;
        }
        var span = _doc == null ? 30 : Math.Max(_doc.Summary().CellA, 10);
        var perPx = span * 1.6 / Math.Max(1, tileHeight) / _orthoZoom;
        _orthoPanX += dx * perPx;
        _orthoPanY -= dy * perPx;
        foreach (var x in ViewportTiles.Where(x => x.Ortho && x.Shown && (_vpLinkZoom || x == t))) RenderViewport(x);
    }

    public void ZoomViewport(ViewportTile t, double wheel)
    {
        var f = Math.Pow(1.12, wheel);
        if (!t.Ortho) { Camera.Zoom = Math.Clamp(Camera.Zoom * f, 0.1, 40); RenderViewport(t); return; }
        _orthoZoom = Math.Clamp(_orthoZoom * f, 0.1, 40);
        foreach (var x in ViewportTiles.Where(x => x.Ortho && x.Shown && (_vpLinkZoom || x == t))) RenderViewport(x);
    }

    /// <summary>A click picks the atom under it in that view (re-rendered at its size first so the pick buffer is its own).</summary>
    public void PickInViewport(ViewportTile t, int px, int py, bool add)
    {
        if (_doc == null) return;
        var opt = ViewportOptions(t);
        var rgba = new byte[t.PixelW * t.PixelH * 4];
        try { _doc.Render(CameraOf(t), opt, rgba); } catch { return; }
        var hit = _doc.Pick(px, py);
        if (hit >= 0) Pick(hit, add); else if (!add) ClearSelection();
        RenderViewports();
    }

    private CapsRenderOpts ViewportOptions(ViewportTile t)
    {
        var o = ViewOptions(t.PixelW, t.PixelH, 1);
        o.Focus = 0;
        return o;
    }

    public void SetViewportSize(ViewportTile t, int w, int h)
    {
        w = Math.Clamp(w, 32, 4096);
        h = Math.Clamp(h, 32, 4096);
        if (w == t.PixelW && h == t.PixelH && t.Image != null) return;
        t.PixelW = w;
        t.PixelH = h;
        RenderViewport(t);
    }

    public void RenderViewports() { foreach (var t in ViewportTiles.Where(t => t.Shown)) RenderViewport(t); }

    private void RenderViewport(ViewportTile t)
    {
        if (_doc == null || !IsViewports) return;
        var doc = _doc;
        var gen = ++t.Gen;
        var cam = CameraOf(t);
        var opt = ViewportOptions(t);
        var (w, h) = (t.PixelW, t.PixelH);
        Task.Run(() =>
        {
            var rgba = new byte[w * h * 4];
            try { doc.Render(cam, opt, rgba); } catch { return; }
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (gen == t.Gen && IsViewports) t.Image = ToBitmap(rgba, w, h); });
        });
    }
}
