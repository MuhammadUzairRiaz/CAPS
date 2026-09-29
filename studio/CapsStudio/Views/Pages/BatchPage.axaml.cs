using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class BatchPage : PageBase
{
    public BatchPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm) return;
            vm.BatchPlotChanged += () =>
            {
                var plot = this.FindControl<LinePlot>("BatchPlot");
                if (plot == null) return;
                plot.YLabel = vm.BatchMetricName;
                var pts = vm.BatchPlotPoints();
                plot.SetData(pts, pts);
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunBatch();
    private void OnStop(object? s, RoutedEventArgs e) => Vm.StopBatch();
    private void OnPause(object? s, RoutedEventArgs e) => Vm.PauseBatch();

    private async void OnAddFiles(object? s, RoutedEventArgs e)
    {
        var top = TopLevel.GetTopLevel(this);
        if (top == null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions { Title = "Inputs for the batch", AllowMultiple = true });
        Vm.AddBatchFiles(files.Select(f => f.TryGetLocalPath()).OfType<string>());
    }

    private async void OnOpenResults(object? s, RoutedEventArgs e)
    {
        if (!File.Exists(Vm.BatchOut)) return;
        var top = TopLevel.GetTopLevel(this);
        if (top != null) await top.Launcher.LaunchFileInfoAsync(new FileInfo(Vm.BatchOut));
    }
}
