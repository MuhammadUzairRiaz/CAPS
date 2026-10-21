using Avalonia.Controls;
using Avalonia.Platform.Storage;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using CapsStudio.ViewModels;

namespace CapsStudio.Views.Pages;

public partial class DiffusionPage : PageBase
{
    private MainViewModel? _hooked;

    public DiffusionPage()
    {
        AvaloniaXamlLoader.Load(this);
        DataContextChanged += (_, _) =>
        {
            if (DataContext is not MainViewModel vm || vm == _hooked) return;
            _hooked = vm;
            vm.DfChanged += () =>
            {
                // the curves are kept as log10 values; drawn on logarithmic axes
                static double P(double v) => Math.Pow(10, v);
                (double, double, string)? window = vm.DfWindow is { } w ? (P(w.From), P(w.To), "fit window") : null;
                var msd = vm.DfMsd;
                this.FindControl<SciPlot>("Msd")!.Figure = msd.Length < 2 ? null : new SciFigure
                {
                    XLabel = "*t* / ps", YLabel = "MSD / Å²", XLog = true, YLog = true, XBand = window,
                    Series = [new SciSeries("chain centres of mass", msd.Select(p => P(p.X)).ToArray(), msd.Select(p => P(p.Y)).ToArray()) { Kind = SciKind.LineMarkers, Width = 1.4 }],
                };
                var sl = vm.DfSlopeCurve;
                this.FindControl<SciPlot>("Slope")!.Figure = sl.Length < 2 ? null : new SciFigure
                {
                    XLabel = "*t* / ps", YLabel = "d ln MSD / d ln *t*", XLog = true, RefY = 1, XBand = window,
                    Series = [new SciSeries("", sl.Select(p => P(p.X)).ToArray(), sl.Select(p => p.Y).ToArray()) { Kind = SciKind.LineMarkers, Width = 1.2 }],
                };
            };
        };
    }

    private async void OnRun(object? s, RoutedEventArgs e) => await Vm.RunDiffusion();
    private void OnCancel(object? s, RoutedEventArgs e) => Vm.Analyze.Cancel();
    private void OnExport(object? s, RoutedEventArgs e) => Window?.ExportAnalysis();
    private void OnRemoveReplica(object? s, RoutedEventArgs e) { if ((s as Control)?.Tag is DfReplica r) Vm.RemoveDfReplica(r); }
    private async void OnAddReplicas(object? s, RoutedEventArgs e)
    {
        if (TopLevel.GetTopLevel(this) is not { } top) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new Avalonia.Platform.Storage.FilePickerOpenOptions
        {
            Title = "Other runs of the same system", AllowMultiple = true,
            FileTypeFilter = [new Avalonia.Platform.Storage.FilePickerFileType("Trajectories") { Patterns = ["*.lammpstrj", "*.dump", "*.dcd", "*.xtc", "*.trr", "*.nc", "*.mdcrd", "*.gro", "*.pdb", "*.xyz"] },
                              new Avalonia.Platform.Storage.FilePickerFileType("All files") { Patterns = ["*"] }],
        });
        Vm.AddDfReplicas(files.Select(f => f.TryGetLocalPath()).OfType<string>());
    }
}
