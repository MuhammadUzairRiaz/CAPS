using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class ComparePage : PageBase
{
    private Point? _drag;

    public ComparePage()
    {
        AvaloniaXamlLoader.Load(this);
        foreach (var (name, a) in new[] { ("TileA", true), ("TileB", false) })
        {
            var tile = this.FindControl<Border>(name)!;
            tile.PointerPressed += (_, e) => { _drag = e.GetPosition(tile); e.Pointer.Capture(tile); };
            tile.PointerMoved += (_, e) =>
            {
                if (_drag is not { } p0) return;
                var p = e.GetPosition(tile);
                Vm.RotateCompare(a, p.X - p0.X, p.Y - p0.Y);
                _drag = p;
            };
            tile.PointerReleased += (_, e) => { _drag = null; e.Pointer.Capture(null); };
        }
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.ComparePlotChanged += () =>
            {
                var plot = this.FindControl<LinePlot>("ComparePlotView");
                if (plot == null) return;
                var (a, b, x, y) = vm.ComparePlot();
                plot.XLabel = x;
                plot.YLabel = y;
                plot.SetCompare(a, b);
            };
        };
    }

    private async Task<string?> Pick(string title)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return null;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = title, AllowMultiple = false });
        return files.Count > 0 ? files[0].TryGetLocalPath() : null;
    }

    private async void OnChooseA(object? s, RoutedEventArgs e) { if (await Pick("Input A") is { } p) await Vm.SetCompareInput(true, p); }
    private async void OnChooseB(object? s, RoutedEventArgs e) { if (await Pick("Input B") is { } p) await Vm.SetCompareInput(false, p); }
    private void OnEditPipeline(object? s, RoutedEventArgs e) => Vm.OpenVisualize();
}
