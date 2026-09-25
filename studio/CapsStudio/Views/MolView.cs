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

    /// <summary>Draws the cell box (a crystal's supercell dashed with one unit cell in the accent).</summary>
    public bool ShowCell { get; set; }
    /// <summary>Drawing style (0 ball and stick, 1 space filling, 2 sticks, 3 no hydrogens, 4 backbone).</summary>
    public int DrawStyle { get; set; }
    /// <summary>Colour mode (0 element, 1 molecule, 2 type, 3 distance to own centre).</summary>
    public int ColourMode { get; set; }
    /// <summary>A dashed line between two atoms with a label at its middle (a chain's end-to-end vector).</summary>
    public int LineA { get; set; } = -1;
    public int LineB { get; set; } = -1;
    public string LineLabel { get; set; } = "";
    /// <summary>Atoms drawn with a selection ring (up to four).</summary>
    public int[] Highlights { get; set; } = [];
    private Point? _pa, _pb;

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
                    Width = w, Height = h, Supersample = scale >= 1.5 ? 1 : 2, Background = light ? 1 : 0, Style = DrawStyle, ColourBy = ColourMode,
                    Outlines = 1, DepthCue = 1, ShowCell = ShowCell ? 1 : 0, Highlight0 = Highlights.Length > 0 ? Highlights[0] : -1, Highlight1 = Highlights.Length > 1 ? Highlights[1] : -1,
                    Highlight2 = Highlights.Length > 2 ? Highlights[2] : -1, Highlight3 = Highlights.Length > 3 ? Highlights[3] : -1,
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
                _pa = _pb = null;
                if (LineA >= 0 && LineB >= 0)
                    try
                    {
                        var n = Math.Max(LineA, LineB) + 1;
                        var p = doc.ProjectAtoms(cam, opt, (int)doc.Summary().Atoms);
                        if (3 * n <= p.Length) { _pa = new Point(p[3 * LineA] / scale, p[3 * LineA + 1] / scale); _pb = new Point(p[3 * LineB] / scale, p[3 * LineB + 1] / scale); }
                    }
                    catch { }
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
        if (_pa is Point a && _pb is Point b)
        {
            var acc = Tokens.Brush("AccB");
            ctx.DrawLine(new Pen(acc, 1.6, new DashStyle([5, 4], 0)), a, b);
            ctx.DrawEllipse(acc, null, a, 3.5, 3.5);
            ctx.DrawEllipse(acc, null, b, 3.5, 3.5);
            if (LineLabel.Length > 0)
            {
                var ft = new FormattedText(LineLabel, System.Globalization.CultureInfo.InvariantCulture, FlowDirection.LeftToRight, new Typeface(Tokens.Mono), 11.5, acc);
                var m = new Point((a.X + b.X) / 2 + 10, (a.Y + b.Y) / 2 - ft.Height / 2);
                ctx.DrawRectangle(new SolidColorBrush(Color.FromArgb(0xE6, 0x16, 0x19, 0x1C)), new Pen(acc, 1), new Rect(m.X - 6, m.Y - 3, ft.Width + 12, ft.Height + 6), 4, 4);
                ctx.DrawText(ft, m);
            }
        }
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
