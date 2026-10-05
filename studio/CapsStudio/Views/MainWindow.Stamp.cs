using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>The clipboard tray, the stamp under the pointer and the two drop targets (design/boards/Stamp).</summary>
public partial class MainWindow
{
    /// <summary>The click's point in the view plane through the structure's centre (the projection fitted to its atoms,
    /// as the move tool does), and the view axis: the stamp goes there, turned about the axis.</summary>
    private void PlaceStampAt(Point at, bool many)
    {
        if (_vm.Document is not { } doc) return;
        var n = (int)doc.Summary().Atoms;
        if (n < 4) { _vm.Status = "The structure needs a few atoms to place a stamp among"; return; }
        var proj = doc.ProjectAtoms(_lastCam, _lastOpt, n);
        var fit = new List<(double X, double Y, double Z, double Sx, double Sy)>();
        var step = Math.Max(1, n / 600);
        for (var i = 0; i < n; i += step)
        {
            if (proj[3 * i + 2] <= 0) continue;
            var a = doc.Atom(i);
            fit.Add((a.X, a.Y, a.Z, proj[3 * i] / _scaling, proj[3 * i + 1] / _scaling));
        }
        if (fit.Count < 4) return;
        // the centre and its place on screen (the fitted map is affine: the mean maps to the mean)
        double cx = fit.Average(f => f.X), cy = fit.Average(f => f.Y), cz = fit.Average(f => f.Z);
        double sx = fit.Average(f => f.Sx), sy = fit.Average(f => f.Sy);
        var d = MainViewModel.ScreenToWorld(fit, at.X - sx, at.Y - sy);
        // the view axis: what a screen move cannot reach (the cross product of two screen-plane directions)
        var u = MainViewModel.ScreenToWorld(fit, 1, 0);
        var v = MainViewModel.ScreenToWorld(fit, 0, 1);
        if (d == null || u == null || v == null) return;
        double[] axis = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        _vm.PlaceStamp([cx + d[0], cy + d[1], cz + d[2]], axis, many);
        RequestRender();
    }

    private void OnTrayPiece(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is TrayPiece p) _vm.ArmStamp(p); ViewHost.Focus(); }
    private void OnTrayRemove(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is TrayPiece p) _vm.RemoveFromTray(p); }
    private void OnTrayHide(object? s, RoutedEventArgs e) { _vm.TrayShown = false; _vm.DisarmStamp(); }

    // ---- Look (design/boards/Look)
    private void OnLookAll(object? s, RoutedEventArgs e) => _vm.LookScope = 0;
    private void OnLookSelection(object? s, RoutedEventArgs e) => _vm.LookScope = 1;
    private void OnLookKeep(object? s, RoutedEventArgs e) => _vm.KeepLook();
    private void OnLookReset(object? s, RoutedEventArgs e) => _vm.ResetLook();

    // ---- dragging a file over the window: two targets, into this structure (a stamp) or open it new
    private void InitDropTargets()
    {
        DragDrop.SetAllowDrop(this, true);
        AddHandler(DragDrop.DragEnterEvent, (_, e) => { if (e.Data.GetFiles() != null && _vm.HasDocument && _vm.IsStudio) _vm.DropTargets = true; });
        AddHandler(DragDrop.DragLeaveEvent, (_, e) =>
        {
            // left the window (not just moved between its children)
            var p = e.GetPosition(this);
            if (p.X <= 0 || p.Y <= 0 || p.X >= Bounds.Width || p.Y >= Bounds.Height) _vm.DropTargets = false;
        });
        AddHandler(DragDrop.DropEvent, (_, _) => _vm.DropTargets = false, Avalonia.Interactivity.RoutingStrategies.Bubble, handledEventsToo: true);
        DragDrop.SetAllowDrop(DropInto, true);
        DragDrop.SetAllowDrop(DropNew, true);
        DropInto.AddHandler(DragDrop.DropEvent, (_, e) =>
        {
            e.Handled = true;
            _vm.DropTargets = false;
            if (e.Data.GetFiles()?.Select(f => f.TryGetLocalPath()).OfType<string>().FirstOrDefault() is { } p) _vm.DropIntoStructure(p);
        });
        DropNew.AddHandler(DragDrop.DropEvent, (_, e) =>
        {
            e.Handled = true;
            _vm.DropTargets = false;
            if (e.Data.GetFiles() is { } files) OpenMany(files.Select(f => f.TryGetLocalPath()).OfType<string>().ToList());
        });
    }
}
