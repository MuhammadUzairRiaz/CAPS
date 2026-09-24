using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Media;
using Avalonia.Media.Imaging;
using Avalonia.Platform;
using CapsStudio.Interop;

namespace CapsStudio.Views;

/// <summary>A small 3D view of a document with its own camera (the molecule builder's preview): ball and stick,
/// coloured by element; drag rotates, the wheel zooms, double-click resets.</summary>
public sealed class MolView : Control
{
    private CapsDocument? _doc;
    private CapsCamera _cam = new() { Yaw = 0.55, Pitch = 0.40, Zoom = 1.0 };
    private WriteableBitmap? _bmp;
    private int _requested, _rendered;
    private bool _busy;
    private Point _last;
    private bool _drag;

    public MolView()
    {
        ClipToBounds = true;
        Focusable = true;
    }

    public CapsDocument? Document
    {
        get => _doc;
        set { _doc = value; Refresh(); }
    }

    public void Reset()
    {
        _cam = new CapsCamera { Yaw = 0.55, Pitch = 0.40, Zoom = 1.0 };
        Refresh();
    }

    public void Refresh()
    {
        _requested++;
        if (!_busy) _ = Loop();
    }

    protected override void OnSizeChanged(SizeChangedEventArgs e)
    {
        base.OnSizeChanged(e);
        Refresh();
    }

    private async Task Loop()
    {
        _busy = true;
        try
        {
            while (_rendered != _requested)
            {
                var ticket = _requested;
                var doc = _doc;
                var scale = VisualRoot?.RenderScaling ?? 1;
                var w = (int)(Math.Max(16, Bounds.Width) * scale);
                var h = (int)(Math.Max(16, Bounds.Height) * scale);
                if (doc == null || Bounds.Width < 16) { _bmp = null; InvalidateVisual(); _rendered = ticket; break; }
                var light = Application.Current?.ActualThemeVariant == Avalonia.Styling.ThemeVariant.Light;
                var opt = new CapsRenderOpts
                {
                    Width = w, Height = h, Supersample = scale >= 1.5 ? 1 : 2, Background = light ? 1 : 0, Style = 0, ColourBy = 0,
                    Outlines = 1, DepthCue = 1, ShowCell = 0, Highlight0 = -1, Highlight1 = -1, Highlight2 = -1, Highlight3 = -1,
                };
                var cam = _cam;
                var buf = new byte[w * h * 4];
                try { await Task.Run(() => doc.Render(cam, opt, buf)); }
                catch { _rendered = ticket; break; }   // disposed while rendering: a newer document follows
                var bmp = new WriteableBitmap(new PixelSize(w, h), new Vector(96 * scale, 96 * scale), PixelFormat.Rgba8888, AlphaFormat.Unpremul);
                using (var fb = bmp.Lock())
                    for (var y = 0; y < h; y++)
                        System.Runtime.InteropServices.Marshal.Copy(buf, y * w * 4, fb.Address + y * fb.RowBytes, w * 4);
                var old = _bmp;
                _bmp = bmp;
                old?.Dispose();
                InvalidateVisual();
                _rendered = ticket;
            }
        }
        finally { _busy = false; }
    }

    public override void Render(DrawingContext ctx)
    {
        ctx.FillRectangle(Tokens.Brush("Bg0B"), new Rect(Bounds.Size));
        if (_bmp != null) ctx.DrawImage(_bmp, new Rect(0, 0, Bounds.Width, Bounds.Height));
    }

    protected override void OnPointerPressed(PointerPressedEventArgs e)
    {
        base.OnPointerPressed(e);
        if (e.ClickCount == 2) { Reset(); return; }
        _last = e.GetPosition(this);
        _drag = true;
        e.Pointer.Capture(this);
    }

    protected override void OnPointerMoved(PointerEventArgs e)
    {
        base.OnPointerMoved(e);
        if (!_drag) return;
        var p = e.GetPosition(this);
        _cam.Yaw += (p.X - _last.X) * 0.01;
        _cam.Pitch = Math.Clamp(_cam.Pitch + (p.Y - _last.Y) * 0.01, -1.5, 1.5);
        _last = p;
        Refresh();
    }

    protected override void OnPointerReleased(PointerReleasedEventArgs e)
    {
        base.OnPointerReleased(e);
        _drag = false;
        e.Pointer.Capture(null);
    }

    protected override void OnPointerWheelChanged(PointerWheelEventArgs e)
    {
        base.OnPointerWheelChanged(e);
        _cam.Zoom = Math.Clamp(_cam.Zoom * Math.Pow(1.12, e.Delta.Y), 0.3, 6);
        Refresh();
    }
}
