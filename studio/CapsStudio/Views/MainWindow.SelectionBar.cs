using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Threading;
using Avalonia.VisualTree;

namespace CapsStudio.Views;

/// <summary>The selection bar over the selection and the ring menu under a held right button (design/boards/SelectionBar).</summary>
public partial class MainWindow
{
    private DispatcherTimer? _barTimer, _ringTimer;
    private bool _ringPending, _ringHold;
    private Point _ringAt;
    private bool _boxMode;
    private (double, double, double, double, double, int) _barCam;   // the camera the bar was placed for

    /// <summary>The Force field page: its table and 3D view take the height left under the header and the controls (at least
    /// 520 px); on a short window the page scrolls instead of squeezing them under the footer.</summary>
    private void FitFieldBody()
    {
        var view = FieldScroll.Viewport.Height;
        if (view <= 0) return;
        var used = FieldHead.Bounds.Height + (FieldControls.IsVisible ? FieldControls.Bounds.Height : 0) + 16 + 12;
        var h = Math.Max(520, Math.Floor(view - used));
        if (Math.Abs(FieldBody.Height - h) > 0.5) FieldBody.Height = h;
    }

    private void InitSelectionBar()
    {
        FieldScroll.PropertyChanged += (_, e) => { if (e.Property == ScrollViewer.ViewportProperty) FitFieldBody(); };
        FieldHead.SizeChanged += (_, _) => FitFieldBody();
        FieldControls.SizeChanged += (_, _) => FitFieldBody();
        _vm.SelectionBarChanged += () => QueuePlaceSelBar(0);
        _vm.PropertyChanged += (_, e) => { if (e.PropertyName is nameof(_vm.HasSelBar)) QueuePlaceSelBar(0); };
        ViewHost.PropertyChanged += (_, e) => { if (e.Property == BoundsProperty) QueuePlaceSelBar(60); };
    }

    /// <summary>The bar is placed once the view is still (never re-projected per frame while it turns).</summary>
    private void QueuePlaceSelBar(int ms)
    {
        if (!_vm.HasSelBar) { SelBar.IsVisible = false; return; }
        _barTimer ??= new DispatcherTimer(TimeSpan.FromMilliseconds(140), DispatcherPriority.Background, (_, _) => { _barTimer!.Stop(); PlaceSelBar(); });
        _barTimer.Stop();
        _barTimer.Interval = TimeSpan.FromMilliseconds(Math.Max(1, ms));
        _barTimer.Start();
    }

    /// <summary>While the view moves the bar steps aside; it returns over the selection when the motion stops.</summary>
    private void SelBarOnViewMoved()
    {
        if (!SelBar.IsVisible && !_vm.HasSelBar) return;
        SelBar.IsVisible = false;
        QueuePlaceSelBar(160);
    }

    private void PlaceSelBar()
    {
        if (!_vm.HasSelBar || _vm.Document is not { } doc || !_haveLast || Ring.IsVisible || _dragging) { SelBar.IsVisible = false; return; }
        float[] xy;
        try
        {
            var sample = _vm.SelBarSample;
            if (_lastFit is { } fit)   // the GPU view: its own camera (the core's last CPU render may be older)
            {
                xy = new float[2 * sample.Length];
                for (var k = 0; k < sample.Length; ++k)
                {
                    var a = doc.Atom(sample[k]);
                    var (px, py, _) = ViewModels.MainViewModel.ProjectFit(fit, (a.X, a.Y, a.Z));
                    xy[2 * k] = (float)px; xy[2 * k + 1] = (float)py;
                }
            }
            else xy = doc.ProjectIndices(_lastCam, _lastOpt, sample);
        }
        catch { SelBar.IsVisible = false; return; }
        var w = ViewHost.Bounds.Width;
        var h = ViewHost.Bounds.Height;
        double x0 = double.MaxValue, y0 = double.MaxValue, x1 = double.MinValue, y1 = double.MinValue;
        var inView = 0;
        for (var k = 0; k + 1 < xy.Length; k += 2)
        {
            var x = xy[k] / _scaling;
            var y = xy[k + 1] / _scaling;
            if (double.IsNaN(x) || double.IsNaN(y) || x < 0 || y < 0 || x > w || y > h) continue;
            ++inView;
            x0 = Math.Min(x0, x); x1 = Math.Max(x1, x); y0 = Math.Min(y0, y); y1 = Math.Max(y1, y);
        }
        if (inView == 0) { SelBar.IsVisible = false; return; }   // the selection is out of the view: the HUD chip still counts it
        SelBar.IsVisible = true;
        SelBar.Measure(Size.Infinity);
        var bw = SelBar.DesiredSize.Width - SelBar.Margin.Left - SelBar.Margin.Right;   // the desired size carries the margin
        var bh = SelBar.DesiredSize.Height - SelBar.Margin.Top - SelBar.Margin.Bottom;
        var cx = (x0 + x1) / 2;
        // above the selection, else below it; always inside the view, clear of the HUD row at the top
        var top = y0 - bh - 16;
        if (top < 44) top = y1 + 16;
        if (top + bh > h - 8) top = Math.Max(44, Math.Min(h - bh - 8, y0 - bh - 16));
        var left = Math.Clamp(cx - bw / 2, 8, Math.Max(8, w - bw - 8));
        SelBar.Margin = new Thickness(left, top, 0, 0);
    }

    private async void OnSelBarAction(object? s, RoutedEventArgs e)
    {
        if ((s as Control)?.Tag is not string id) return;
        if (id == "more") { OpenSelBarMore(s as Control); return; }
        await RunSelAction(id);
    }

    private async Task RunSelAction(string id)
    {
        var text = await _vm.RunSelectionAction(id);
        if (text is { Length: > 0 } && Clipboard != null) await Clipboard.SetTextAsync(text);
        ViewHost.Focus();
    }

    private void OpenSelBarMore(Control? anchor)
    {
        var items = new List<MenuItem>();
        MenuItem M(string header, string id, string gesture = "")
        {
            var m = new MenuItem { Header = header, InputGesture = gesture.Length > 0 ? KeyGesture.Parse(gesture) : null };
            m.Click += async (_, _) => await RunSelAction(id);
            return m;
        }
        items.Add(M("Frame the selection", "frame", "F"));
        items.Add(M("Show only the selection", "only"));
        items.Add(M("Invert the selection", "invert"));
        var grow = new MenuItem { Header = "Grow along bonds" };
        grow.Click += (_, _) => { _vm.GrowSelectionKey(); ViewHost.Focus(); };
        items.Add(grow);
        if (_vm.HasHiddenAtoms) items.Add(M("Show all atoms", "showall"));
        var menu = new ContextMenu { ItemsSource = items };
        menu.Open(anchor ?? SelBar);
    }

    private void OnShowAllAtoms(object? s, RoutedEventArgs e) { _vm.ShowAllAtoms(); ViewHost.Focus(); }

    // ---------------------------------------------------------------- layers (design/boards/Layers)

    private void OnLayerRow(object? s, TappedEventArgs e)
    {
        if ((s as Control)?.Tag is not ViewModels.LayerRow r || e.Source is Visual v && v.FindAncestorOfType<Button>() != null) return;
        _vm.SelectLayer(r, e.KeyModifiers.HasFlag(KeyModifiers.Shift));
    }
    private void OnLayerExpand(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.LayerRow r) _vm.ToggleLayer(r); e.Handled = true; }
    private void OnLayerEye(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.LayerRow r) _vm.CycleLayer(r); e.Handled = true; }
    private void OnLayerLock(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is ViewModels.LayerRow r) _vm.LockLayer(r); e.Handled = true; }

    // ---------------------------------------------------------------- the ring

    /// <summary>A right press: after a short hold without moving, the ring opens under the pointer (a drag pans as before).</summary>
    private void RingPressed(Point at)
    {
        _ringPending = true;
        _ringAt = at;
        _ringTimer ??= new DispatcherTimer(TimeSpan.FromMilliseconds(280), DispatcherPriority.Input, (_, _) =>
        {
            _ringTimer!.Stop();
            if (_ringPending && _dragging && !_moved) OpenRing(_ringAt, hold: true);
        });
        _ringTimer.Stop();
        _ringTimer.Start();
    }

    private void OpenRing(Point at, bool hold)
    {
        _ringPending = false;
        _ringTimer?.Stop();
        if (_vm.Document is not { } doc || _vm.Busy || !_vm.IsStudio) return;
        // nothing selected and an atom under the pointer: the ring acts on it
        if (!_vm.HasSelBar)
        {
            var hit = PickAtView(doc, at);
            if (hit >= 0) { _vm.StartSelection(hit, "click"); _vm.RefreshSelBar(); RequestRender(); }
        }
        _ringHold = hold;
        SelBar.IsVisible = false;
        // wholly inside the view (near an edge the ring moves in; aiming is by direction, so the flick still works)
        var o = RingMenu.Outer + 2;
        var c = new Point(Math.Clamp(at.X, o, Math.Max(o, ViewHost.Bounds.Width - o)), Math.Clamp(at.Y, o, Math.Max(o, ViewHost.Bounds.Height - o)));
        Ring.Open(c, _vm.RingItems(), _vm.HasSelBar ? _vm.SelBarCountText : "");
        _vm.Status = hold ? "Point at an action and let go (anywhere in the middle: nothing)" : "Click an action (Esc closes)";
    }

    /// <summary>For the screenshot tool: the ring a held right button opens.</summary>
    public void HoldRingForTest(Point window) { if (!Ring.IsVisible && _dragging && !_moved && this.TranslatePoint(window, ViewHost) is { } at) OpenRing(at, hold: true); }

    private void CloseRing()
    {
        Ring.Close();
        _ringPending = false;
        QueuePlaceSelBar(0);
    }

    /// <summary>The ring's pointer move: aims at a wedge. True when the ring owns the pointer.</summary>
    private bool RingMoved(Point at)
    {
        if (_ringPending && (Math.Abs(at.X - _ringAt.X) + Math.Abs(at.Y - _ringAt.Y) > 4)) _ringPending = false;   // a drag: pan
        if (!Ring.IsVisible) return false;
        Ring.Aim(at);
        return true;
    }

    /// <summary>The release that ends a held ring runs the aimed action; a quick right click leaves the ring open.</summary>
    private bool RingReleased(Point at, bool rightClickWithoutMove)
    {
        if (_ringPending && rightClickWithoutMove) { OpenRing(at, hold: false); return true; }
        _ringPending = false;
        if (!Ring.IsVisible) return false;
        if (!_ringHold) return true;   // the sticky ring waits for a click
        var hot = Ring.Aim(at);
        var items = Ring.Items;
        CloseRing();
        if (hot >= 0 && hot < items.Count) _ = RunSelAction(items[hot].Id);
        else _vm.Status = "";
        return true;
    }

    /// <summary>A click while the sticky ring is open: its action, or (outside it, or in the hub) nothing.</summary>
    private bool RingClick(Point at)
    {
        if (!Ring.IsVisible) return false;
        var d = at - Ring.Centre;
        var r = Math.Sqrt(d.X * d.X + d.Y * d.Y);
        var hot = r <= RingMenu.Outer ? Ring.Aim(at) : -1;
        var items = Ring.Items;
        CloseRing();
        if (hot >= 0 && hot < items.Count) _ = RunSelAction(items[hot].Id);
        return true;
    }

    // ---------------------------------------------------------------- the box

    /// <summary>The box as the lasso's polygon (its four corners), so the same inside test selects.</summary>
    private void UpdateBox(Point at)
    {
        _lassoPts = [_press, new Point(at.X, _press.Y), at, new Point(_press.X, at.Y)];
        Labels.SetLasso(_lassoPts);
    }
}
